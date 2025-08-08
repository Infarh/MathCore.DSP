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
        
        // Вычисляем массив фаз
        var phases = new float[samples.Length];
        for (var i = 0; i < samples.Length; i++)
        {
            phases[i] = GetPhase(samples[i]);
        }
        
        // Разворачиваем фазы (phase unwrapping)
        UnwrapPhases(phases);
        
        // Вычисляем мгновенную частоту как производную фазы
        var dt = 1.0 / fd;
        var scale_factor = (float)(fd / (2.0 * Math.PI));
        
        for (var i = 1; i < samples.Length; i++)
        {
            // Производная фазы дает мгновенную частоту в рад/с
            // Делим на 2π для перевода в Гц
            var instantaneous_frequency = (phases[i] - phases[i - 1]) * scale_factor;
            
            // Вычитаем центральную частоту f0
            result[i] = instantaneous_frequency - (float)f0;
        }
        
        return result;
    }

    /// <summary>Быстрое вычисление фазы с использованием статических функций SampleSI16</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static float GetPhase(SampleSI16 sample)
    {
        // Используем уже оптимизированную функцию GetArg из SampleSI16
        return SampleSI16.GetArg(sample);
    }

    /// <summary>Разворачивание фаз - устранение скачков ±2π</summary>
    private static void UnwrapPhases(Span<float> phases)
    {
        if (phases.Length <= 1) return;
        
        const float two_pi = (float)(2.0 * Math.PI);
        const float pi = (float)Math.PI;
        
        for (var i = 1; i < phases.Length; i++)
        {
            var diff = phases[i] - phases[i - 1];
            
            // Приводим разность к диапазону (-π, π]
            while (diff > pi)
            {
                diff -= two_pi;
                phases[i] -= two_pi;
            }
            while (diff <= -pi)
            {
                diff += two_pi;
                phases[i] += two_pi;
            }
        }
    }
}
