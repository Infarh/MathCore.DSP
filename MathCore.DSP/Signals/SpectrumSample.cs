namespace MathCore.DSP.Signals;

/// <summary>Отсчёт сигнала</summary>
/// <remarks>Новый отсчёт сигнала</remarks>
/// <param name="Frequency">Частота</param>
/// <param name="Value">Значение</param>
public readonly struct SpectrumSample(double Frequency, Complex Value) : IEquatable<SpectrumSample>
{
    /// <summary>Частота</summary>
    public double Frequency { get; init; } = Frequency;

    /// <summary>Значение</summary>
    public Complex Value { get; init; } = Value;

    public bool Equals(SpectrumSample other) => Frequency == other.Frequency && Value == other.Value;

    public override bool Equals(object? obj) => obj is SignalSample other && Equals(other);

    public override int GetHashCode() => unchecked((Frequency.GetHashCode() * 397) ^ Value.GetHashCode());

    public override string ToString() => $"{Frequency}:{Value}";

    public static bool operator ==(SpectrumSample left, SpectrumSample right) => left.Equals(right);

    public static bool operator !=(SpectrumSample left, SpectrumSample right) => !left.Equals(right);

    public static implicit operator Complex(SpectrumSample sample) => sample.Value;

    public void Deconstruct(out double frequency, out Complex value) => (frequency, value) = (Frequency, Value);
}