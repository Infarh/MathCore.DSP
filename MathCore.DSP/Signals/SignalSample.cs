// ReSharper disable UnusedType.Global
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Signals;

/// <summary>Отсчёт сигнала</summary>
/// <remarks>Новый отсчёт сигнала</remarks>
/// <param name="Time">Время</param>
/// <param name="Value">Значение</param>
public readonly struct SignalSample(double Time, double Value) : IEquatable<SignalSample>
{
    /// <summary>Новый отсчёт сигнала</summary>
    public SignalSample() : this(0, 0) { }

    /// <summary>Время</summary>
    public double Time { get; init; } = Time;

    /// <summary>Значение</summary>
    public double Value { get; init; } = Value;

    public bool Equals(SignalSample other) => Time == other.Time && Value == other.Value;

    public override bool Equals(object? obj) => obj is SignalSample other && Equals(other);

    public override int GetHashCode() => unchecked((Time.GetHashCode() * 397) ^ Value.GetHashCode());

    public override string ToString() => $"{Time}:{Value}";

    public static bool operator ==(SignalSample left, SignalSample right) => left.Equals(right);

    public static bool operator !=(SignalSample left, SignalSample right) => !left.Equals(right);

    public static implicit operator double(SignalSample sample) => sample.Value;

    public void Deconstruct(out double Time, out double Value) => (Time, Value) = (this.Time, this.Value);
}