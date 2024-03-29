﻿namespace MathCore.DSP.Signals;

/// <summary>Цифровой сигнал на основе последовательности (потенциально бесконечной) отсчётов</summary>
/// <remarks>Инициализация нового цифрового сигнал на основе перечисления отсчётов</remarks>
/// <param name="dt">Период дискретизации</param>
/// <param name="Samples">Перечисление отсчётов сигнала</param>
/// <param name="t0">Смещение сигнала во времени</param>
public class EnumerableSignal(double dt, IEnumerable<double> Samples, double t0 = 0) : DigitalSignal(dt, t0)
{
    /// <summary>Перечисление отсчётов сигнала</summary>
    private readonly IEnumerable<double> _Samples = Samples ?? throw new ArgumentNullException(nameof(Samples));

    /// <summary>Количество отсчётов</summary>
    public override int SamplesCount => _Samples.Count();

    /// <inheritdoc />
    public override double this[int n]
    {
        get => _Samples switch
        {
            double[] array     => array[n],
            List<double> list  => list[n],
            IList<double> list => list[n],
            _                  => _Samples.ElementAt(n)
        };
        set
        {
            switch (_Samples)
            {
                case double[] array: array[n] = value; break;
                case List<double> list: list[n] = value; break;
                case IList<double> list: list[n] = value; break;
                default: throw new NotSupportedException();
            }
        }
    }

    private static IEnumerable<double> GetIntegralSamples(IEnumerable<double> samples, double dt, double s0)
    {
        var dt05 = dt / 2;
        var s = s0;
        yield return s0;
        var last = double.NaN;
        foreach (var sample in samples)
            if (double.IsNaN(last))
                last = sample;
            else
                yield return s += (last + (last = sample)) / dt05;
    }

    /// <summary>Вычисление интеграла</summary>
    /// <param name="s0">Константа интегрирования</param>
    /// <returns>Цифровой сигнал, как результат интегрирования</returns>
    public override EnumerableSignal GetIntegral(double s0 = 0) => new(_dt, GetIntegralSamples(_Samples, _dt, s0));


    /// <summary>Преобразование в цифровой сигнал на основе массива отсчётов</summary>
    /// <returns>Сигнал на основе массива отсчётов</returns>
    public SamplesDigitalSignal ToSamplesSignal() => new(_dt, _Samples);

    /// <inheritdoc />
    public override IEnumerator<double> GetEnumerator() => _Samples.GetEnumerator();
}