// ReSharper disable UnusedType.Global
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Signals;

public readonly struct SignalSample : IEquatable<SignalSample>
{
    public double Time { get; }

    public double Value { get; }

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

public static class SignalSampleEx
{
    public static IEnumerable<double> AsDouble(this IEnumerable<SignalSample> samples) => samples.Select(s => s.Value);

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