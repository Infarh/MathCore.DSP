
using MathCore.HackRF;
using MathCore.HackRF.Streaming;

const int capture_seconds = 2;
const double sample_rate_hz = 10_000_000;
const uint center_frequency_hz = 90_000_000;
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
