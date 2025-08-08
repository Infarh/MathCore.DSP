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
}
