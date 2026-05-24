#nullable enable

using ComplexN = System.Numerics.Complex;

namespace MathCore.DSP.SDR;

/// <summary>Набор инструментов цифровой обработки SDR-сигналов для IQ и аудио</summary>
public static class SdrSignalProcessing
{
    /// <summary>Рассчитать базовые метрики interleaved IQ-потока int8</summary>
    /// <param name="RawIq">Байтовый буфер формата I,Q,I,Q</param>
    /// <returns>Кортеж средних, RMS, доли клиппинга и средней мощности</returns>
    /// <exception cref="ArgumentNullException">Буфер исходных данных не задан</exception>
    /// <exception cref="ArgumentException">Длина буфера не кратна двум</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var metrics = SdrSignalProcessing.ComputeIq8Metrics(raw_iq);
    /// Console.WriteLine($"Mean I={metrics.mean_i:F4}, Mean Q={metrics.mean_q:F4}");
    /// ]]></code>
    /// </remarks>
    public static (double mean_i, double mean_q, double rms_i, double rms_q, double clip_ratio, double avg_power) ComputeIq8Metrics(byte[] RawIq)
    {
        ArgumentNullException.ThrowIfNull(RawIq);

        if ((RawIq.Length & 1) != 0)
            throw new ArgumentException("Размер массива IQ должен быть кратен двум", nameof(RawIq));

        if (RawIq.Length == 0)
            return (0, 0, 0, 0, 0, 0);

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

    /// <summary>Рассчитать средние значения I и Q и мощность DC-компоненты комплексного сигнала</summary>
    /// <param name="Samples">Массив комплексных отсчётов</param>
    /// <returns>Кортеж средних I/Q и мощности DC</returns>
    /// <exception cref="ArgumentNullException">Массив отсчётов не задан</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var dc_metrics = SdrSignalProcessing.ComputeComplexDcMetrics(iq_after_dc);
    /// Console.WriteLine($"Pdc={dc_metrics.dc_power:E6}");
    /// ]]></code>
    /// </remarks>
    public static (double mean_i, double mean_q, double dc_power) ComputeComplexDcMetrics(ComplexN[] Samples)
    {
        ArgumentNullException.ThrowIfNull(Samples);

        if (Samples.Length == 0)
            return (0, 0, 0);

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

    /// <summary>Преобразовать отношение мощностей в децибелы</summary>
    /// <param name="Ratio">Отношение мощностей</param>
    /// <returns>Значение в dB</returns>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var dc_suppression_db = SdrSignalProcessing.ToDb(stage2_before_pdc / stage2_after_pdc);
    /// ]]></code>
    /// </remarks>
    public static double ToDb(double Ratio) => 10 * Math.Log10(Math.Max(Ratio, 1e-18));

    /// <summary>Убрать DC-компоненту из interleaved IQ-потока int8</summary>
    /// <param name="RawIq">Байтовый буфер формата I,Q,I,Q</param>
    /// <param name="Alpha">Коэффициент сглаживания экспоненциальной оценки среднего в диапазоне (0;1]</param>
    /// <returns>Массив комплексных IQ-отсчётов после подавления DC</returns>
    /// <exception cref="ArgumentNullException">Буфер исходных данных не задан</exception>
    /// <exception cref="ArgumentOutOfRangeException">Коэффициент сглаживания выходит за допустимый диапазон</exception>
    /// <exception cref="ArgumentException">Длина буфера не кратна двум</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var iq_after_dc = SdrSignalProcessing.RemoveDcFromInterleavedIq8(raw_iq, Alpha: 0.0025);
    /// ]]></code>
    /// </remarks>
    public static ComplexN[] RemoveDcFromInterleavedIq8(byte[] RawIq, double Alpha)
    {
        ArgumentNullException.ThrowIfNull(RawIq);

        if (Alpha <= 0 || Alpha > 1)
            throw new ArgumentOutOfRangeException(nameof(Alpha), Alpha, "Коэффициент сглаживания должен быть в диапазоне (0;1]");

        if ((RawIq.Length & 1) != 0)
            throw new ArgumentException("Размер массива IQ должен быть кратен двум", nameof(RawIq));

        var count = RawIq.Length / 2;
        var result = new ComplexN[count];

        var mean_i = 0d;
        var mean_q = 0d;

        for (var index = 0; index < count; index++)
        {
            var i_value = unchecked((sbyte)RawIq[2 * index]);
            var q_value = unchecked((sbyte)RawIq[2 * index + 1]);

            var i_d = (double)i_value;
            var q_d = (double)q_value;

            mean_i += Alpha * (i_d - mean_i);
            mean_q += Alpha * (q_d - mean_q);

            result[index] = new(i_d - mean_i, q_d - mean_q);
        }

        return result;
    }

    /// <summary>Выполнить частотный сдвиг комплексного IQ-сигнала цифровым гетеродином</summary>
    /// <param name="Source">Исходный комплексный сигнал</param>
    /// <param name="ShiftHz">Частота сдвига в герцах</param>
    /// <param name="SampleRateHz">Частота дискретизации входного сигнала в герцах</param>
    /// <returns>Новый массив IQ-отсчётов после сдвига частоты</returns>
    /// <exception cref="ArgumentNullException">Исходный массив не задан</exception>
    /// <exception cref="ArgumentOutOfRangeException">Частота дискретизации не является положительной</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var shifted_iq = SdrSignalProcessing.ShiftFrequency(iq_after_dc, ShiftHz: 800_000, SampleRateHz: 10_000_000);
    /// ]]></code>
    /// </remarks>
    public static ComplexN[] ShiftFrequency(ComplexN[] Source, double ShiftHz, double SampleRateHz)
    {
        ArgumentNullException.ThrowIfNull(Source);

        if (SampleRateHz <= 0)
            throw new ArgumentOutOfRangeException(nameof(SampleRateHz), SampleRateHz, "Частота дискретизации должна быть положительной");

        var result = new ComplexN[Source.Length];
        var dphase = -2 * Math.PI * ShiftHz / SampleRateHz;

        for (var index = 0; index < Source.Length; index++)
        {
            var phase = dphase * index;
            var lo = ComplexN.FromPolarCoordinates(1, phase);
            result[index] = Source[index] * lo;
        }

        return result;
    }

    /// <summary>Оценить мощность узкого канала вокруг заданной частоты через комплексное смешение и однополюсный НЧ-фильтр</summary>
    /// <param name="Source">Входной IQ-сигнал</param>
    /// <param name="SampleRateHz">Частота дискретизации входного сигнала</param>
    /// <param name="CenterHz">Центральная частота канала для оценки мощности</param>
    /// <param name="BandwidthHz">Эффективная полоса оценочного НЧ-фильтра</param>
    /// <returns>Средняя мощность канала</returns>
    /// <exception cref="ArgumentOutOfRangeException">Некорректные параметры частоты дискретизации или полосы</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var p_center = SdrSignalProcessing.EstimateChannelPower(iq_data, SampleRateHz: 10_000_000, CenterHz: 0, BandwidthHz: 150_000);
    /// var p_offset = SdrSignalProcessing.EstimateChannelPower(iq_data, SampleRateHz: 10_000_000, CenterHz: 800_000, BandwidthHz: 150_000);
    /// ]]></code>
    /// </remarks>
    public static double EstimateChannelPower(ReadOnlySpan<ComplexN> Source, double SampleRateHz, double CenterHz, double BandwidthHz)
    {
        if (SampleRateHz <= 0)
            throw new ArgumentOutOfRangeException(nameof(SampleRateHz), SampleRateHz, "Частота дискретизации должна быть положительной");

        if (BandwidthHz <= 0)
            throw new ArgumentOutOfRangeException(nameof(BandwidthHz), BandwidthHz, "Полоса оценочного фильтра должна быть положительной");

        if (Source.Length == 0)
            return 0;

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

    /// <summary>Синтезировать КИХ-ФНЧ оконным методом с Blackman-окном</summary>
    /// <param name="SampleRateHz">Частота дискретизации сигнала в герцах</param>
    /// <param name="PassbandHz">Частота конца полосы пропускания в герцах</param>
    /// <param name="TransitionHz">Ширина переходной полосы в герцах</param>
    /// <param name="TapCount">Количество коэффициентов фильтра</param>
    /// <returns>Массив коэффициентов КИХ-фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Некорректные параметры синтеза фильтра</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var taps = SdrSignalProcessing.DesignLowPassFirBlackman(
    ///     SampleRateHz: 10_000_000,
    ///     PassbandHz: 75_000,
    ///     TransitionHz: 25_000,
    ///     TapCount: 401);
    /// ]]></code>
    /// </remarks>
    public static double[] DesignLowPassFirBlackman(double SampleRateHz, double PassbandHz, double TransitionHz, int TapCount)
    {
        if (SampleRateHz <= 0)
            throw new ArgumentOutOfRangeException(nameof(SampleRateHz), SampleRateHz, "Частота дискретизации должна быть положительной");

        if (PassbandHz <= 0)
            throw new ArgumentOutOfRangeException(nameof(PassbandHz), PassbandHz, "Частота полосы пропускания должна быть положительной");

        if (TransitionHz <= 0)
            throw new ArgumentOutOfRangeException(nameof(TransitionHz), TransitionHz, "Переходная полоса должна быть положительной");

        if (TapCount <= 3 || (TapCount & 1) == 0)
            throw new ArgumentOutOfRangeException(nameof(TapCount), TapCount, "Количество коэффициентов должно быть нечётным и больше трёх");

        var taps = new double[TapCount];
        var middle = TapCount / 2;
        var cutoff_hz = PassbandHz + TransitionHz * 0.5;
        var fc = cutoff_hz / SampleRateHz;

        for (var n = 0; n < TapCount; n++)
        {
            var k = n - middle;
            var sinc = k == 0
                ? 2 * fc
                : Math.Sin(2 * Math.PI * fc * k) / (Math.PI * k);

            var window = 0.42
                       - 0.5 * Math.Cos(2 * Math.PI * n / (TapCount - 1))
                       + 0.08 * Math.Cos(4 * Math.PI * n / (TapCount - 1));

            taps[n] = sinc * window;
        }

        var sum = 0d;
        for (var n = 0; n < taps.Length; n++)
            sum += taps[n];

        if (Math.Abs(sum) > 1e-18)
            for (var n = 0; n < taps.Length; n++)
                taps[n] /= sum;

        return taps;
    }

    /// <summary>Выполнить КИХ-фильтрацию комплексного сигнала и децимацию</summary>
    /// <param name="Source">Входной IQ-сигнал</param>
    /// <param name="Taps">Коэффициенты КИХ-фильтра</param>
    /// <param name="Decimation">Коэффициент децимации</param>
    /// <returns>Отфильтрованный и децимированный сигнал</returns>
    /// <exception cref="ArgumentNullException">Входные массивы не заданы</exception>
    /// <exception cref="ArgumentOutOfRangeException">Коэффициент децимации должен быть положительным</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var iq_400k = SdrSignalProcessing.LowPassDecimateComplex(shifted_iq, taps, Decimation: 25);
    /// ]]></code>
    /// </remarks>
    public static ComplexN[] LowPassDecimateComplex(ComplexN[] Source, double[] Taps, int Decimation)
    {
        ArgumentNullException.ThrowIfNull(Source);
        ArgumentNullException.ThrowIfNull(Taps);

        if (Decimation <= 0)
            throw new ArgumentOutOfRangeException(nameof(Decimation), Decimation, "Коэффициент децимации должен быть положительным");

        if (Source.Length == 0 || Taps.Length == 0)
            return [];

        var half = Taps.Length / 2;
        var output = new List<ComplexN>(Source.Length / Decimation + 1);

        for (var n = half; n < Source.Length - half; n += Decimation)
        {
            var i_acc = 0d;
            var q_acc = 0d;

            for (var k = -half; k <= half; k++)
            {
                var tap = Taps[k + half];
                var sample = Source[n - k];
                i_acc += sample.Real * tap;
                q_acc += sample.Imaginary * tap;
            }

            output.Add(new(i_acc, q_acc));
        }

        return output.ToArray();
    }

    /// <summary>Выполнить рациональный ресемплинг комплексного сигнала по windowed-sinc ядру</summary>
    /// <param name="Source">Входной IQ-сигнал</param>
    /// <param name="Interpolation">Коэффициент интерполяции</param>
    /// <param name="Decimation">Коэффициент децимации</param>
    /// <param name="HalfTaps">Половина длины sinc-ядра</param>
    /// <returns>Сигнал после пересчёта частоты дискретизации</returns>
    /// <exception cref="ArgumentNullException">Входной массив не задан</exception>
    /// <exception cref="ArgumentOutOfRangeException">Коэффициенты пересчёта заданы некорректно</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var iq_960k = SdrSignalProcessing.ResampleRationalSinc(iq_400k, Interpolation: 12, Decimation: 5, HalfTaps: 16);
    /// ]]></code>
    /// </remarks>
    public static ComplexN[] ResampleRationalSinc(ComplexN[] Source, int Interpolation, int Decimation, int HalfTaps)
    {
        ArgumentNullException.ThrowIfNull(Source);

        if (Interpolation <= 0)
            throw new ArgumentOutOfRangeException(nameof(Interpolation), Interpolation, "Коэффициент интерполяции должен быть положительным");

        if (Decimation <= 0)
            throw new ArgumentOutOfRangeException(nameof(Decimation), Decimation, "Коэффициент децимации должен быть положительным");

        if (HalfTaps < 2)
            throw new ArgumentOutOfRangeException(nameof(HalfTaps), HalfTaps, "Размер половины sinc-ядра должен быть не меньше двух");

        if (Source.Length == 0)
            return [];

        var output_length = (int)Math.Floor(Source.Length * (double)Interpolation / Decimation);
        var output = new ComplexN[output_length];

        for (var out_index = 0; out_index < output_length; out_index++)
        {
            var src_position = out_index * (double)Decimation / Interpolation;
            var src_center = (int)Math.Floor(src_position);
            var frac = src_position - src_center;

            var i_acc = 0d;
            var q_acc = 0d;
            var w_sum = 0d;

            for (var tap = -HalfTaps; tap <= HalfTaps; tap++)
            {
                var src_index = src_center + tap;
                if (src_index < 0 || src_index >= Source.Length)
                    continue;

                var x = tap - frac;
                var sinc = Math.Abs(x) < 1e-12
                    ? 1d
                    : Math.Sin(Math.PI * x) / (Math.PI * x);

                var window = 0.54 + 0.46 * Math.Cos(Math.PI * x / (HalfTaps + 1));
                var weight = sinc * window;

                var sample = Source[src_index];
                i_acc += sample.Real * weight;
                q_acc += sample.Imaginary * weight;
                w_sum += weight;
            }

            if (Math.Abs(w_sum) > 1e-15)
            {
                i_acc /= w_sum;
                q_acc /= w_sum;
            }

            output[out_index] = new(i_acc, q_acc);
        }

        return output;
    }

    /// <summary>Выполнить FM-демодуляцию и децимацию аудиосигнала</summary>
    /// <param name="Source">Входной комплексный сигнал</param>
    /// <param name="SourceSampleRateHz">Частота дискретизации входного сигнала</param>
    /// <param name="Decimation">Коэффициент децимации на выходе демодулятора</param>
    /// <param name="PassbandHz">Частота полосы пропускания аудио-ФНЧ</param>
    /// <param name="TransitionHz">Переходная полоса аудио-ФНЧ</param>
    /// <param name="TapCount">Количество коэффициентов аудио-ФНЧ</param>
    /// <returns>Массив отсчётов демодулированного аудиосигнала</returns>
    /// <exception cref="ArgumentNullException">Входной массив не задан</exception>
    /// <exception cref="ArgumentOutOfRangeException">Параметры фильтрации или децимации заданы некорректно</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var audio_96k = SdrSignalProcessing.FmDemodulateAndDecimate(
    ///     iq_960k,
    ///     SourceSampleRateHz: 960_000,
    ///     Decimation: 10,
    ///     PassbandHz: 45_000,
    ///     TransitionHz: 10_000,
    ///     TapCount: 161);
    /// ]]></code>
    /// </remarks>
    public static double[] FmDemodulateAndDecimate(
        ComplexN[] Source,
        double SourceSampleRateHz,
        int Decimation,
        double PassbandHz,
        double TransitionHz,
        int TapCount)
    {
        ArgumentNullException.ThrowIfNull(Source);

        if (SourceSampleRateHz <= 0)
            throw new ArgumentOutOfRangeException(nameof(SourceSampleRateHz), SourceSampleRateHz, "Частота дискретизации должна быть положительной");

        if (Decimation <= 0)
            throw new ArgumentOutOfRangeException(nameof(Decimation), Decimation, "Коэффициент децимации должен быть положительным");

        if (Source.Length < 2)
            return [];

        var demodulated = new double[Source.Length - 1];
        var scale = SourceSampleRateHz / (2 * Math.PI);

        for (var index = 1; index < Source.Length; index++)
        {
            var cross = ComplexN.Conjugate(Source[index - 1]) * Source[index];
            var phase = Math.Atan2(cross.Imaginary, cross.Real);
            demodulated[index - 1] = phase * scale;
        }

        var taps = DesignLowPassFirBlackman(SourceSampleRateHz, PassbandHz, TransitionHz, TapCount);
        return LowPassDecimateReal(demodulated, taps, Decimation);
    }

    /// <summary>Выполнить КИХ-фильтрацию вещественного сигнала и децимацию</summary>
    /// <param name="Source">Входной массив вещественных отсчётов</param>
    /// <param name="Taps">Коэффициенты КИХ-фильтра</param>
    /// <param name="Decimation">Коэффициент децимации</param>
    /// <returns>Отфильтрованный и децимированный вещественный сигнал</returns>
    /// <exception cref="ArgumentNullException">Входные массивы не заданы</exception>
    /// <exception cref="ArgumentOutOfRangeException">Коэффициент децимации должен быть положительным</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var audio_48k = SdrSignalProcessing.LowPassDecimateReal(audio_96k, audio_lpf_taps, Decimation: 2);
    /// ]]></code>
    /// </remarks>
    public static double[] LowPassDecimateReal(double[] Source, double[] Taps, int Decimation)
    {
        ArgumentNullException.ThrowIfNull(Source);
        ArgumentNullException.ThrowIfNull(Taps);

        if (Decimation <= 0)
            throw new ArgumentOutOfRangeException(nameof(Decimation), Decimation, "Коэффициент децимации должен быть положительным");

        if (Source.Length == 0 || Taps.Length == 0)
            return [];

        var half = Taps.Length / 2;
        var output = new List<double>(Source.Length / Decimation + 1);

        for (var n = half; n < Source.Length - half; n += Decimation)
        {
            var acc = 0d;

            for (var k = -half; k <= half; k++)
                acc += Source[n - k] * Taps[k + half];

            output.Add(acc);
        }

        return output.ToArray();
    }

    /// <summary>Рассчитать RMS вещественного сигнала</summary>
    /// <param name="Source">Входной массив отсчётов</param>
    /// <returns>Среднеквадратичное значение сигнала</returns>
    /// <exception cref="ArgumentNullException">Входной массив не задан</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var rms_value = SdrSignalProcessing.ComputeRms(audio_data);
    /// ]]></code>
    /// </remarks>
    public static double ComputeRms(double[] Source)
    {
        ArgumentNullException.ThrowIfNull(Source);

        if (Source.Length == 0)
            return 0;

        var sum = 0d;
        for (var i = 0; i < Source.Length; i++)
            sum += Source[i] * Source[i];

        return Math.Sqrt(sum / Source.Length);
    }

    /// <summary>Нормализовать аудиосигнал по пику</summary>
    /// <param name="Source">Входной вещественный сигнал</param>
    /// <param name="TargetPeak">Целевой пик амплитуды после нормализации</param>
    /// <returns>Массив float-отсчётов после нормализации</returns>
    /// <exception cref="ArgumentNullException">Входной массив не задан</exception>
    /// <exception cref="ArgumentOutOfRangeException">Целевой пик должен быть положительным</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var normalized_audio = SdrSignalProcessing.NormalizePeak(audio_data, TargetPeak: 0.8);
    /// ]]></code>
    /// </remarks>
    public static float[] NormalizePeak(double[] Source, double TargetPeak)
    {
        ArgumentNullException.ThrowIfNull(Source);

        if (TargetPeak <= 0)
            throw new ArgumentOutOfRangeException(nameof(TargetPeak), TargetPeak, "Целевой пик должен быть положительным");

        var max_abs = 0d;
        for (var i = 0; i < Source.Length; i++)
            max_abs = Math.Max(max_abs, Math.Abs(Source[i]));

        var scale = max_abs > 1e-12 ? TargetPeak / max_abs : 1;
        var result = new float[Source.Length];

        for (var i = 0; i < Source.Length; i++)
            result[i] = (float)(Source[i] * scale);

        return result;
    }
}
