
using MathCore.HackRF;
using MathCore.HackRF.Streaming;
using System.Numerics;

const int capture_seconds = 2;
const double sample_rate_hz = 10_000_000;
const uint center_frequency_hz = 90_000_000;
const double station_offset_hz = 800_000;
const uint lna_gain_db = 24;
const uint vga_gain_db = 24;

Console.WriteLine("Этап 1. Захват IQ с HackRF");

Device.Initialize();

try
{
	var captured_data = CaptureIqBlock(
		capture_seconds,
		sample_rate_hz,
		center_frequency_hz,
		lna_gain_db,
		vga_gain_db,
		out var rx_statistics);

	var stage1_metrics = ComputeIqMetrics(captured_data);

	Console.WriteLine($"Получено байт: {captured_data.Length}");
	Console.WriteLine($"Получено IQ-отсчётов: {captured_data.Length / 2}");
	Console.WriteLine($"Mean I: {stage1_metrics.mean_i:F4}");
	Console.WriteLine($"Mean Q: {stage1_metrics.mean_q:F4}");
	Console.WriteLine($"RMS I: {stage1_metrics.rms_i:F2}");
	Console.WriteLine($"RMS Q: {stage1_metrics.rms_q:F2}");
	Console.WriteLine($"Clip ratio: {stage1_metrics.clip_ratio * 100:F3}%");
	Console.WriteLine($"Power avg: {stage1_metrics.avg_power:F2}");
	Console.WriteLine($"Rx dropped blocks: {rx_statistics.DroppedBlocks}");
	Console.WriteLine($"Rx max processing ms: {rx_statistics.MaxProcessingMilliseconds:F3}");
	Console.WriteLine($"Этап 1 завершён");

	Console.WriteLine();
	Console.WriteLine("Этап 2. Подавление DC-пика");

	var iq_after_dc = RemoveDcPeak(captured_data, alpha: 0.0025);

	var stage2_before = ComputeComplexMetrics(iq_after_dc.before_dc);
	var stage2_after = ComputeComplexMetrics(iq_after_dc.after_dc);

	Console.WriteLine($"DC до: I={stage2_before.mean_i:F5}, Q={stage2_before.mean_q:F5}, Pdc={stage2_before.dc_power:F6}");
	Console.WriteLine($"DC после: I={stage2_after.mean_i:F5}, Q={stage2_after.mean_q:F5}, Pdc={stage2_after.dc_power:F6}");
	Console.WriteLine($"Ослабление DC (dB): {ToDb(stage2_before.dc_power / Math.Max(stage2_after.dc_power, 1e-18)):F2}");
	Console.WriteLine("Этап 2 завершён");

	Console.WriteLine();
	Console.WriteLine("Этап 3. Гетеродин +0.8 МГц");

	var shifted_iq = MixByHeterodyne(iq_after_dc.after_dc, station_offset_hz, sample_rate_hz);

	var validation_count = Math.Min(1_000_000, shifted_iq.Length);
	var before_slice = iq_after_dc.after_dc.AsSpan(0, validation_count);
	var after_slice = shifted_iq.AsSpan(0, validation_count);

	var channel_power_before = EstimateChannelPower(before_slice, sample_rate_hz, station_offset_hz, 150_000);
	var zero_power_before = EstimateChannelPower(before_slice, sample_rate_hz, 0, 150_000);

	var channel_power_after = EstimateChannelPower(after_slice, sample_rate_hz, station_offset_hz, 150_000);
	var zero_power_after = EstimateChannelPower(after_slice, sample_rate_hz, 0, 150_000);

	Console.WriteLine($"До сдвига: P(800кГц)={channel_power_before:E6}, P(0)={zero_power_before:E6}");
	Console.WriteLine($"После сдвига: P(800кГц)={channel_power_after:E6}, P(0)={zero_power_after:E6}");
	Console.WriteLine($"Усиление целевого канала в нуле (dB): {ToDb(zero_power_after / Math.Max(zero_power_before, 1e-18)):F2}");
	Console.WriteLine($"Подавление 800кГц после сдвига (dB): {ToDb(channel_power_before / Math.Max(channel_power_after, 1e-18)):F2}");
	Console.WriteLine("Этап 3 завершён");
}
finally
{
	Device.Shutdown();
}

Console.WriteLine("End.");
return;

static byte[] CaptureIqBlock(
	int CaptureSeconds,
	double SampleRateHz,
	uint CenterFrequencyHz,
	uint LnaGainDb,
	uint VgaGainDb,
	out RxSessionStatistics Statistics)
{
	using var device = new Device();

	device.SampleRate = SampleRateHz;
	device.Frequency = CenterFrequencyHz;
	device.LnaGain = LnaGainDb;
	device.VgaGain = VgaGainDb;
	device.EnableLNA = true;
	device.FilterBandwidth = 10_000_000;

	var target_iq_samples = (int)(CaptureSeconds * SampleRateHz);
	var target_bytes = target_iq_samples * 2;
	var result = new byte[target_bytes];
	var received_bytes = 0;

	using var data_ready = new ManualResetEventSlim(false);

	using var rx_session = device.StartRxSession(
		(ReadOnlySpan<byte> rx_block, in RxBlockMetadata _) =>
		{
			if (received_bytes >= target_bytes)
				return;

			var remain_bytes = target_bytes - received_bytes;
			var copy_length = Math.Min(remain_bytes, rx_block.Length);
			rx_block[..copy_length].CopyTo(result.AsSpan(received_bytes, copy_length));
			received_bytes += copy_length;

			if (received_bytes >= target_bytes)
				data_ready.Set();
		},
		new RxSessionOptions
		{
			QueueCapacity = 128,
			OverflowPolicy = RxQueueOverflowPolicy.DropOldest,
			OnProcessingError = error => Console.WriteLine($"RX processing error: {error.Message}")
		});

	if (!data_ready.Wait(TimeSpan.FromSeconds(CaptureSeconds + 5)))
		throw new TimeoutException("Не удалось получить требуемый объём IQ-данных за отведённое время");

	Statistics = rx_session.GetStatistics();
	return result;
}

static (double mean_i, double mean_q, double rms_i, double rms_q, double clip_ratio, double avg_power) ComputeIqMetrics(byte[] RawIq)
{
	var i_sum = 0d;
	var q_sum = 0d;
	var i_pow_sum = 0d;
	var q_pow_sum = 0d;
	var pwr_sum = 0d;
	var clip_count = 0;

	var count = RawIq.Length / 2;
	for (var index = 0; index < RawIq.Length; index += 2)
	{
		var i_value = unchecked((sbyte)RawIq[index]);
		var q_value = unchecked((sbyte)RawIq[index + 1]);

		i_sum += i_value;
		q_sum += q_value;

		i_pow_sum += i_value * i_value;
		q_pow_sum += q_value * q_value;
		pwr_sum += i_value * i_value + q_value * q_value;

		if (i_value is sbyte.MinValue or sbyte.MaxValue)
			clip_count++;
		if (q_value is sbyte.MinValue or sbyte.MaxValue)
			clip_count++;
	}

	return (
		mean_i: i_sum / count,
		mean_q: q_sum / count,
		rms_i: Math.Sqrt(i_pow_sum / count),
		rms_q: Math.Sqrt(q_pow_sum / count),
		clip_ratio: (double)clip_count / (count * 2),
		avg_power: pwr_sum / count);
}

static (Complex[] before_dc, Complex[] after_dc) RemoveDcPeak(byte[] RawIq, double alpha)
{
	var count = RawIq.Length / 2;
	var before_dc = new Complex[count];
	var after_dc = new Complex[count];

	var mean_i = 0d;
	var mean_q = 0d;

	for (var i = 0; i < count; i++)
	{
		var i_value = unchecked((sbyte)RawIq[2 * i]);
		var q_value = unchecked((sbyte)RawIq[2 * i + 1]);

		var i_d = (double)i_value;
		var q_d = (double)q_value;

		before_dc[i] = new Complex(i_d, q_d);

		mean_i += alpha * (i_d - mean_i);
		mean_q += alpha * (q_d - mean_q);

		after_dc[i] = new Complex(i_d - mean_i, q_d - mean_q);
	}

	return (before_dc, after_dc);
}

static (double mean_i, double mean_q, double dc_power) ComputeComplexMetrics(Complex[] Samples)
{
	var i_sum = 0d;
	var q_sum = 0d;

	for (var i = 0; i < Samples.Length; i++)
	{
		i_sum += Samples[i].Real;
		q_sum += Samples[i].Imaginary;
	}

	var mean_i = i_sum / Samples.Length;
	var mean_q = q_sum / Samples.Length;
	return (mean_i, mean_q, mean_i * mean_i + mean_q * mean_q);
}

static double ToDb(double Ratio) => 10 * Math.Log10(Math.Max(Ratio, 1e-18));

static Complex[] MixByHeterodyne(Complex[] Source, double ShiftHz, double SampleRateHz)
{
	var result = new Complex[Source.Length];
	var dphase = -2 * Math.PI * ShiftHz / SampleRateHz;

	for (var index = 0; index < Source.Length; index++)
	{
		var phase = dphase * index;
		var lo = Complex.FromPolarCoordinates(1, phase);
		result[index] = Source[index] * lo;
	}

	return result;
}

static double EstimateChannelPower(ReadOnlySpan<Complex> Source, double SampleRateHz, double CenterHz, double BandwidthHz)
{
	var phase_step = -2 * Math.PI * CenterHz / SampleRateHz;
	var alpha = Math.Min(1, 2 * Math.PI * BandwidthHz / SampleRateHz);

	var i_lp = 0d;
	var q_lp = 0d;
	var pwr_sum = 0d;

	for (var index = 0; index < Source.Length; index++)
	{
		var phase = phase_step * index;
		var lo_i = Math.Cos(phase);
		var lo_q = Math.Sin(phase);

		var sample_i = Source[index].Real;
		var sample_q = Source[index].Imaginary;

		var mix_i = sample_i * lo_i - sample_q * lo_q;
		var mix_q = sample_i * lo_q + sample_q * lo_i;

		i_lp += alpha * (mix_i - i_lp);
		q_lp += alpha * (mix_q - q_lp);

		pwr_sum += i_lp * i_lp + q_lp * q_lp;
	}

	return pwr_sum / Source.Length;
}
