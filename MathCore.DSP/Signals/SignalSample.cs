using System;

// ReSharper disable UnusedType.Global

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Signals
{
    public readonly struct SignalSample : IEquatable<SignalSample>
    {
        public readonly double Time;
        public readonly double Value;

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

        public static bool operator ==(SignalSample left, SignalSample right) => left.Equals(right);
        public static bool operator !=(SignalSample left, SignalSample right) => !left.Equals(right);
    }
}
