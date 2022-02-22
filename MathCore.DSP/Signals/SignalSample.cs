// ReSharper disable UnusedType.Global
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Signals;

/// <summary>Отсчёт сигнала</summary>
public readonly struct SignalSample : IEquatable<SignalSample>
{
    /// <summary>Время</summary>
    public double Time { get; init; }

    /// <summary>Значение</summary>
    public double Value { get; init; }

    /// <summary>Новый отсчёт сигнала</summary>
    /// <param name="Time">Время</param>
    /// <param name="Value">Значение</param>
    public SignalSample(double Time, double Value)
    {
        this.Time = Time;
        this.Value = Value;
    }

    public bool Equals(SignalSample other) => Time.Equals(other.Time) && Value.Equals(other.Value);

    public override bool Equals(object obj) => obj is SignalSample other && Equals(other);

    public override int GetHashCode()
    {
        unchecked
        {
            return (Time.GetHashCode() * 397) ^ Value.GetHashCode();
        }
    }

    public override string ToString() => $"{Time}:{Value}";

    public static bool operator ==(SignalSample left, SignalSample right) => left.Equals(right);
    public static bool operator !=(SignalSample left, SignalSample right) => !left.Equals(right);

    public static implicit operator double(SignalSample sample) => sample.Value;
}

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