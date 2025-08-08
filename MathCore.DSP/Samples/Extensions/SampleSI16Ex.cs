using System.Runtime.CompilerServices;

namespace MathCore.DSP.Samples.Extensions;

public static class SampleSI16Ex
{
    /// <summary>Фазовая демодуляция радиосигнала</summary>
    /// <param name="samples">Последовательность отсчётов квадратурного радиосигнала</param>
    /// <param name="f0">Центральная частота фазовой модуляции</param>
    /// <param name="fd">Частота дискретизации</param>
    /// <returns>Возвращает массив вещественных значений отсчётов демодулированного сигнала</returns>
    public static float[] PhaseDemodulation(this Span<SampleSI16> samples, double f0, double fd)
    {
        if (samples.IsEmpty) return [];
        if (samples.Length == 1) return [0f];

        var result = new float[samples.Length];
        result[0] = 0f; // Первый отсчёт всегда 0, так как нет предыдущего для вычисления производной

        var scale_factor = (float)(fd / (2.0 * Math.PI));
        var f0_offset = (float)f0;
        var unwrapped_correction = 0f;

        for (var i = 1; i < samples.Length; i++)
        {
            // Вычисляем разность фаз напрямую через комплексное произведение
            // Δφ = arg(z_curr * z_prev*) = atan2(Im(z_curr * z_prev*), Re(z_curr * z_prev*))
            var curr = samples[i];
            var prev = samples[i - 1];

            // z_curr * z_prev* = (I1 + jQ1) * (I2 - jQ2) = (I1*I2 + Q1*Q2) + j(Q1*I2 - I1*Q2)
            var real_part = curr.I * prev.I + curr.Q * prev.Q;
            var imag_part = curr.Q * prev.I - curr.I * prev.Q;

            // Разность фаз с одним вызовом atan2
#if NET8_0_OR_GREATER
            var phase_diff = MathF.Atan2(imag_part, real_part);
#else
            var phase_diff = (float)Math.Atan2(imag_part, real_part);
#endif

            // Unwrapping - приводим к диапазону (-π, π] и накапливаем коррекцию
            phase_diff = UnwrapSinglePhaseDiff(phase_diff, ref unwrapped_correction);

            // Мгновенная частота = производная фазы / (2π)
            var instantaneous_frequency = phase_diff * scale_factor;

            // Вычитаем центральную частоту
            result[i] = instantaneous_frequency - f0_offset;
        }

        return result;
    }

    /// <summary>Фазовая модуляция сигнала для передачи</summary>
    /// <param name="data">Модулирующие данные (мгновенные частоты в Гц)</param>
    /// <param name="f0">Центральная частота передачи</param>
    /// <param name="fd">Частота дискретизации</param>
    /// <param name="amplitude">Амплитуда выходного сигнала (по умолчанию 120 для максимального использования динамического диапазона)</param>
    /// <returns>Массив квадратурных отсчётов для передачи</returns>
    public static SampleSI16[] PhaseModulation(this ReadOnlySpan<float> data, double f0, double fd, float amplitude = 120f)
    {
        if (data.IsEmpty) return [];

        var result = new SampleSI16[data.Length];
        var dt = 1.0 / fd;
        var omega0 = 2.0 * Math.PI * f0;
        var accumulated_phase = 0.0;

        for (var i = 0; i < data.Length; i++)
        {
            // Мгновенная частота = f0 + данные[i] 
            var instantaneous_frequency = f0 + data[i];

            // Интегрируем частоту для получения фазы: φ(t) = ∫ω(t)dt
            var omega_instant = 2.0 * Math.PI * instantaneous_frequency;
            accumulated_phase += omega_instant * dt;

            // Генерируем квадратурный сигнал: I = A*cos(φ), Q = A*sin(φ)
#if NET8_0_OR_GREATER
            var cos_phase = MathF.Cos((float)accumulated_phase);
            var sin_phase = MathF.Sin((float)accumulated_phase);
#else
            var cos_phase = (float)Math.Cos(accumulated_phase);
            var sin_phase = (float)Math.Sin(accumulated_phase);
#endif

            // Масштабируем и ограничиваем в диапазоне sbyte
            var i_sample = ClampToSByte(amplitude * cos_phase);
            var q_sample = ClampToSByte(amplitude * sin_phase);

            result[i] = new SampleSI16(i_sample, q_sample);
        }

        return result;
    }

    /// <summary>Фазовая модуляция с опорным сигналом для непрерывности фазы</summary>
    /// <param name="data">Модулирующие данные (мгновенные частоты в Гц)</param>
    /// <param name="f0">Центральная частота передачи</param>
    /// <param name="fd">Частота дискретизации</param>
    /// <param name="initial_phase">Начальная фаза для обеспечения непрерывности</param>
    /// <param name="amplitude">Амплитуда выходного сигнала</param>
    /// <returns>Массив квадратурных отсчётов и финальная фаза для следующего блока</returns>
    public static (SampleSI16[] samples, double final_phase) PhaseModulation(
        this ReadOnlySpan<float> data,
        double f0,
        double fd,
        double initial_phase,
        float amplitude = 120f)
    {
        if (data.IsEmpty) return ([], initial_phase);

        var result = new SampleSI16[data.Length];
        var dt = 1.0 / fd;
        var accumulated_phase = initial_phase;

        for (var i = 0; i < data.Length; i++)
        {
            // Мгновенная частота = f0 + данные[i] 
            var instantaneous_frequency = f0 + data[i];

            // Интегрируем частоту для получения фазы
            var omega_instant = 2.0 * Math.PI * instantaneous_frequency;
            accumulated_phase += omega_instant * dt;

            // Генерируем квадратурный сигнал
#if NET8_0_OR_GREATER
            var cos_phase = MathF.Cos((float)accumulated_phase);
            var sin_phase = MathF.Sin((float)accumulated_phase);
#else
            var cos_phase = (float)Math.Cos(accumulated_phase);
            var sin_phase = (float)Math.Sin(accumulated_phase);
#endif

            var i_sample = ClampToSByte(amplitude * cos_phase);
            var q_sample = ClampToSByte(amplitude * sin_phase);

            result[i] = new SampleSI16(i_sample, q_sample);
        }

        return (result, accumulated_phase);
    }

    /// <summary>Unwrapping одной разности фаз с накоплением коррекции</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static float UnwrapSinglePhaseDiff(float phase_diff, ref float accumulated_correction)
    {
#if NET8_0_OR_GREATER
        const float pi = MathF.PI;
#else
        const float pi = (float)Math.PI;
#endif
        const float two_pi = 2f * pi;

        // Приводим разность к диапазону (-π, π]
        if (phase_diff > pi)
        {
            accumulated_correction -= two_pi;
            phase_diff -= two_pi;
        }
        else if (phase_diff <= -pi)
        {
            accumulated_correction += two_pi;
            phase_diff += two_pi;
        }

        return phase_diff;
    }

    /// <summary>Ограничивает значение в диапазоне sbyte с оптимальным использованием динамического диапазона</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static sbyte ClampToSByte(float value) =>
        value switch // Ограничиваем в диапазоне [-127, 127] для избежания переполнения
        {
            > 127f => 127,
            < -128f => -128,
            _ => (sbyte)Math.Round(value)
        };
}
