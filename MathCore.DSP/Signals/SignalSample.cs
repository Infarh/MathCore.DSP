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
    public SignalSample() { }

    /// <summary>Новый отсчёт сигнала</summary>
    /// <param name="Time">Время</param>
    /// <param name="Value">Значение</param>
    public SignalSample(double Time, double Value) => (this.Time, this.Value) = (Time, Value);

    public bool Equals(SignalSample other) => Time == other.Time && Value == other.Value;

    public override bool Equals(object obj) => obj is SignalSample other && Equals(other);

    public override int GetHashCode() => unchecked((Time.GetHashCode() * 397) ^ Value.GetHashCode());

    public override string ToString() => $"{Time}:{Value}";

    public static bool operator ==(SignalSample left, SignalSample right) => left.Equals(right);

    public static bool operator !=(SignalSample left, SignalSample right) => !left.Equals(right);

    public static implicit operator double(SignalSample sample) => sample.Value;

    public void Deconstruct(out double Time, out double Value) => (Time, Value) = (this.Time, this.Value);
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

/// <summary>Отсчёт сигнала</summary>
public readonly struct SpectrumSample : IEquatable<SpectrumSample>
{
    /// <summary>Частота</summary>
    public double Frequency { get; init; }

    /// <summary>Значение</summary>
    public Complex Value { get; init; }

    /// <summary>Новый отсчёт сигнала</summary>
    /// <param name="Frequency">Частота</param>
    /// <param name="Value">Значение</param>
    public SpectrumSample(double Frequency, Complex Value)
    {
        this.Frequency = Frequency;
        this.Value = Value;
    }

    public bool Equals(SpectrumSample other) => Frequency == other.Frequency && Value == other.Value;

    public override bool Equals(object obj) => obj is SignalSample other && Equals(other);

    public override int GetHashCode() => unchecked((Frequency.GetHashCode() * 397) ^ Value.GetHashCode());

    public override string ToString() => $"{Frequency}:{Value}";

    public static bool operator ==(SpectrumSample left, SpectrumSample right) => left.Equals(right);

    public static bool operator !=(SpectrumSample left, SpectrumSample right) => !left.Equals(right);

    public static implicit operator Complex(SpectrumSample sample) => sample.Value;

    public void Deconstruct(out double Frequency, out Complex Value) => (Frequency, Value) = (this.Frequency, this.Value);
}


/// <summary>Класс методов-расширений для отсчётов спектра</summary>
public static class SpectrumSampleEx
{
    /// <summary>Представление последовательности отсчётов в виде последовательности комплексных чисел - значений</summary>
    /// <param name="samples">Последовательность отсчётов спектра</param>
    /// <returns>Последовательность комплексных чисел - отсчётов спектра</returns>
    public static IEnumerable<Complex> AsComplex(this IEnumerable<SpectrumSample> samples) => samples.Select(s => s.Value);

    /// <summary>Представление последовательности комплексных чисел в виде последовательности отсчётов спектра</summary>
    /// <param name="samples">Последовательность комплексных чисел - значений отсчётов спектра</param>
    /// <param name="df">Разрешающая способность спектра</param>
    /// <param name="f0">Начальное значение частоты - смещение</param>
    /// <returns>Последовательность отсчётов спектра</returns>
    public static IEnumerable<SpectrumSample> AsSamples(this IEnumerable<Complex> samples, double df, double f0 = 0)
    {
        var index = 0;
        foreach (var sample in samples)
        {
            yield return new(index * df + f0, sample);
            index++;
        }
    }
}