﻿namespace MathCore.DSP.Signals;

/// <summary>Класс методов-расширений для отсчётов сигнала</summary>
public static class SignalSampleEx
{
    /// <summary>Представление последовательности отсчётов в виде последовательности вещественных чисел - значений</summary>
    /// <param name="samples">Последовательность отсчётов сигнала</param>
    /// <returns>Последовательность вещественных чисел - отсчётов сигнала</returns>
    public static IEnumerable<double> AsDouble(this IEnumerable<SignalSample> samples) => samples.Select(s => s.Value);

    /// <summary>Представление последовательности вещественных чисел в виде последовательности отсчётов сигнала</summary>
    /// <param name="samples">Последовательность вещественных чисел - значений отсчётов сигнала</param>
    /// <param name="dt">Период времени между отсчётами - период дискретизации</param>
    /// <param name="t0">Начальный момент времени - смещение</param>
    /// <returns>Последовательность отсчётов сигнала</returns>
    public static IEnumerable<SignalSample> AsSamples(this IEnumerable<double> samples, double dt, double t0 = 0)
    {
        var index = 0;
        foreach (var sample in samples)
        {
            yield return new(index * dt + t0, sample);
            index++;
        }
    }
}