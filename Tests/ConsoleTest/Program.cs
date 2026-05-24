
using MathCore.HackRF;
using MathCore.HackRF.Streaming;
using MathCore.DSP.SDR;
using NAudio.Wave;
using System.Diagnostics;
using System.Numerics;

const int capture_seconds = 1;
const double sample_rate_hz = 10_000_000;
const uint center_frequency_hz = 90_000_000;
const double station_offset_hz = 800_000;
const double stage4_passband_hz = 75_000;
const double stage4_transition_hz = 25_000;
const int stage4_decimation = 25;
const int stage5_interpolation = 12;
const int stage5_decimation = 5;
const int stage6_decimation = 10;
const int stage7_decimation = 2;
const uint lna_gain_db = 24;
const uint vga_gain_db = 24;

Console.WriteLine("Этап 1. Захват IQ с HackRF");

Device.Initialize(); // Явно инициализируем native-библиотеку HackRF перед работой

try
{
    // Читаем фиксированный по длительности блок IQ для воспроизводимой диагностики
    var captured_data = CaptureIqBlock(
        capture_seconds,
        sample_rate_hz,
        center_frequency_hz,
        lna_gain_db,
        vga_gain_db,
        out var rx_statistics);

    var stage1_metrics = SdrSignalProcessing.ComputeIQ8Metrics(captured_data); // Базовые метрики качества входного тракта

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

    var iq_after_dc = SdrSignalProcessing.RemoveDcFromInterleavedIQ8(captured_data, Alpha: 0.0025); // Подавляем пик на нулевой частоте

    var stage2_before = (mean_i: stage1_metrics.mean_i, mean_q: stage1_metrics.mean_q, dc_power: stage1_metrics.mean_i * stage1_metrics.mean_i + stage1_metrics.mean_q * stage1_metrics.mean_q);
    var stage2_after = SdrSignalProcessing.ComputeComplexDCMetrics(iq_after_dc);

    Console.WriteLine($"DC до: I={stage2_before.mean_i:F5}, Q={stage2_before.mean_q:F5}, Pdc={stage2_before.dc_power:F6}");
    Console.WriteLine($"DC после: I={stage2_after.mean_i:F5}, Q={stage2_after.mean_q:F5}, Pdc={stage2_after.dc_power:F6}");
    Console.WriteLine($"Ослабление DC (dB): {SdrSignalProcessing.ToDb(stage2_before.dc_power / Math.Max(stage2_after.dc_power, 1e-18)):F2}");
    Console.WriteLine("Этап 2 завершён");

    Console.WriteLine();
    Console.WriteLine("Этап 3. Гетеродин +0.8 МГц");

    var shifted_iq = SdrSignalProcessing.ShiftFrequency(iq_after_dc, ShiftHz: station_offset_hz, SampleRateHz: sample_rate_hz); // Переносим станцию 90.8 МГц в ноль

    var validation_count = Math.Min(1_000_000, shifted_iq.Length); // Берём ограниченный фрагмент для быстрой валидации мощности канала
    var before_slice = iq_after_dc.AsSpan(0, validation_count);
    var after_slice = shifted_iq.AsSpan(0, validation_count);

    var channel_power_before = SdrSignalProcessing.EstimateChannelPower(before_slice, sample_rate_hz, station_offset_hz, 150_000);
    var zero_power_before = SdrSignalProcessing.EstimateChannelPower(before_slice, sample_rate_hz, 0, 150_000);

    var channel_power_after = SdrSignalProcessing.EstimateChannelPower(after_slice, sample_rate_hz, station_offset_hz, 150_000);
    var zero_power_after = SdrSignalProcessing.EstimateChannelPower(after_slice, sample_rate_hz, 0, 150_000);

#if DEBUG
    // В debug-режиме дополнительно проверяем перенос канала через FFT-пики
    var (freq_before_hz, power_before_db) = SdrDebugSpectrum.GetSpectrumSnapshot(iq_after_dc, sample_rate_hz, FftSize: 32768);
    var (freq_after_hz, power_after_db) = SdrDebugSpectrum.GetSpectrumSnapshot(shifted_iq, sample_rate_hz, FftSize: 32768);
    var peak_before = SdrDebugSpectrum.FindPeakInRange(freq_before_hz, power_before_db, 600_000, 1_000_000);
    var peak_after = SdrDebugSpectrum.FindPeakInRange(freq_after_hz, power_after_db, -100_000, 100_000);

    Console.WriteLine($"DEBUG FFT peak до гетеродина: f={peak_before.frequency_hz:F0} Гц, p={peak_before.power_db:F2} dB");
    Console.WriteLine($"DEBUG FFT peak после гетеродина: f={peak_after.frequency_hz:F0} Гц, p={peak_after.power_db:F2} dB");
#endif

    Console.WriteLine($"До сдвига: P(800кГц)={channel_power_before:E6}, P(0)={zero_power_before:E6}");
    Console.WriteLine($"После сдвига: P(800кГц)={channel_power_after:E6}, P(0)={zero_power_after:E6}");
    Console.WriteLine($"Усиление целевого канала в нуле (dB): {SdrSignalProcessing.ToDb(zero_power_after / Math.Max(zero_power_before, 1e-18)):F2}");
    Console.WriteLine($"Подавление 800кГц после сдвига (dB): {SdrSignalProcessing.ToDb(channel_power_before / Math.Max(channel_power_after, 1e-18)):F2}");
    Console.WriteLine("Этап 3 завершён");

    Console.WriteLine();
    Console.WriteLine("Этап 4. ФНЧ + децимация x25");

    var stage4_sample_rate_hz = sample_rate_hz / stage4_decimation;
    var stage4_taps = SdrSignalProcessing.DesignLowPassFirBlackman(sample_rate_hz, stage4_passband_hz, stage4_transition_hz, TapCount: 401); // Антиалиасинговый ФНЧ перед decimate x25
    var stage4_iq = SdrSignalProcessing.LowPassDecimateComplex(shifted_iq, stage4_taps, stage4_decimation);

    var stage4_center_power = SdrSignalProcessing.EstimateChannelPower(stage4_iq, stage4_sample_rate_hz, 0, 75_000);
    var stage4_edge_power = SdrSignalProcessing.EstimateChannelPower(stage4_iq, stage4_sample_rate_hz, 170_000, 20_000);

    Console.WriteLine($"Выходная частота дискретизации: {stage4_sample_rate_hz:F0} Гц");
    Console.WriteLine($"Выходных IQ-отсчётов: {stage4_iq.Length}");
    Console.WriteLine($"P(center)= {stage4_center_power:E6}");
    Console.WriteLine($"P(edge 170кГц)= {stage4_edge_power:E6}");
    Console.WriteLine($"Center/Edge (dB): {SdrSignalProcessing.ToDb(stage4_center_power / Math.Max(stage4_edge_power, 1e-18)):F2}");
    Console.WriteLine("Этап 4 завершён");

    Console.WriteLine();
    Console.WriteLine("Этап 5. Ресемплинг x12/x5");

    var stage5_sample_rate_hz = stage4_sample_rate_hz * stage5_interpolation / stage5_decimation;
    var stage5_iq = SdrSignalProcessing.ResampleRationalSinc(stage4_iq, stage5_interpolation, stage5_decimation, HalfTaps: 16); // Рациональный ресемплинг 12/5

    var stage5_center_power = SdrSignalProcessing.EstimateChannelPower(stage5_iq, stage5_sample_rate_hz, 0, 75_000);
    var stage5_edge_power = SdrSignalProcessing.EstimateChannelPower(stage5_iq, stage5_sample_rate_hz, 300_000, 30_000);

    Console.WriteLine($"Выходная частота дискретизации: {stage5_sample_rate_hz:F0} Гц");
    Console.WriteLine($"Выходных IQ-отсчётов: {stage5_iq.Length}");
    Console.WriteLine($"P(center)= {stage5_center_power:E6}");
    Console.WriteLine($"P(edge 300кГц)= {stage5_edge_power:E6}");
    Console.WriteLine($"Center/Edge (dB): {SdrSignalProcessing.ToDb(stage5_center_power / Math.Max(stage5_edge_power, 1e-18)):F2}");
    Console.WriteLine("Этап 5 завершён");

    Console.WriteLine();
    Console.WriteLine("Этап 6. FM-демодуляция + децимация x10");

    var stage6_sample_rate_hz = stage5_sample_rate_hz / stage6_decimation;
    var stage6_audio = SdrSignalProcessing.FmDemodulateAndDecimate(stage5_iq, stage5_sample_rate_hz, stage6_decimation, 45_000, 10_000, 161); // Фазовый дискриминатор + ФНЧ аудио

    var stage6_rms = SdrSignalProcessing.ComputeRms(stage6_audio);
    Console.WriteLine($"Частота дискретизации после демодуляции: {stage6_sample_rate_hz:F0} Гц");
    Console.WriteLine($"Отсчётов аудио: {stage6_audio.Length}");
    Console.WriteLine($"RMS аудио (96к): {stage6_rms:F2}");
    Console.WriteLine("Этап 6 завершён");

    Console.WriteLine();
    Console.WriteLine("Этап 7. Подготовка 48 кГц и вывод на динамики");

    var stage7_sample_rate_hz = stage6_sample_rate_hz / stage7_decimation;
    var stage7_taps = SdrSignalProcessing.DesignLowPassFirBlackman(stage6_sample_rate_hz, 15_000, 5_000, 127); // Ограничиваем аудиополосу перед финальной децимацией
    var stage7_audio = SdrSignalProcessing.LowPassDecimateReal(stage6_audio, stage7_taps, stage7_decimation);

    var stage7_audio_norm = SdrSignalProcessing.NormalizePeak(stage7_audio, 0.8);
    Console.WriteLine($"Частота дискретизации финального аудио: {stage7_sample_rate_hz:F0} Гц");
    Console.WriteLine($"Отсчётов финального аудио: {stage7_audio_norm.Length}");

    PlayAudio(stage7_audio_norm, (int)stage7_sample_rate_hz);
    Console.WriteLine("Этап 7 завершён");
}
finally
{
    Device.Shutdown(); // Гарантированно освобождаем ресурсы HackRF даже при исключениях
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

    // Параметры радиотракта для сценария приёма FM-станции
    device.SampleRate = SampleRateHz;
    device.Frequency = CenterFrequencyHz;
    device.LnaGain = LnaGainDb;
    device.VgaGain = VgaGainDb;
    device.EnableLNA = true;
    device.FilterBandwidth = 10_000_000;

    var target_iq_samples = (int)(CaptureSeconds * SampleRateHz); // Количество комплексных отсчётов для целевой длительности
    var target_bytes = target_iq_samples * 2;
    var result = new byte[target_bytes];
    var received_bytes = 0;

    using var data_ready = new ManualResetEventSlim(false); // Сигнал о наборе полного блока данных

    using var rx_session = device.StartRxSession(
        (ReadOnlySpan<byte> rx_block, in RxBlockMetadata _) =>
        {
            if (received_bytes >= target_bytes)
                return;

            // Копируем только нужный остаток, чтобы не выйти за пределы целевого буфера
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
            OverflowPolicy = RxQueueOverflowPolicy.DropOldest, // Предпочитаем свежие данные при перегрузке
            OnProcessingError = error => Console.WriteLine($"RX processing error: {error.Message}")
        });

    if (!data_ready.Wait(TimeSpan.FromSeconds(CaptureSeconds + 5)))
        throw new TimeoutException("Не удалось получить требуемый объём IQ-данных за отведённое время");

    Statistics = rx_session.GetStatistics();
    return result;
}

static void PlayAudio(float[] Audio, int SampleRate)
{
    // Формируем одноканальный PCM16-поток для NAudio
    var wave_format = new WaveFormat(SampleRate, 16, 1);
    using var wave_out = new WaveOutEvent();
    var provider = new BufferedWaveProvider(wave_format)
    {
        DiscardOnBufferOverflow = true
    };

    var pcm16 = new byte[Audio.Length * 2];
    for (var i = 0; i < Audio.Length; i++)
    {
        // Нормализованные float-отсчёты [-1;1] преобразуем в signed 16-bit little-endian
        var sample = (short)Math.Clamp(Math.Round(Audio[i] * short.MaxValue), short.MinValue, short.MaxValue);
        pcm16[2 * i] = (byte)(sample & 0xff);
        pcm16[2 * i + 1] = (byte)((sample >> 8) & 0xff);
    }

    provider.AddSamples(pcm16, 0, pcm16.Length);
    wave_out.Init(provider);
    wave_out.Play();

    var playback_timeout_ms = (int)(Audio.Length * 1000d / SampleRate) + 3000; // Таймаут защищает от зависания устройства вывода
    var timer = Stopwatch.StartNew();

    while ((wave_out.PlaybackState == PlaybackState.Playing || provider.BufferedBytes > 0) && timer.ElapsedMilliseconds < playback_timeout_ms)
        Thread.Sleep(50);

    wave_out.Stop();
}
