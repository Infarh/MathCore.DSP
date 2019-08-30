using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Signals
{
    public class SamplesDigitalSignal : DigitalSignal, IEquatable<SamplesDigitalSignal>
    {
        [NotNull] private readonly double[] _Samples;

        public override int SamplesCount => _Samples.Length;

        public override double this[int n]
        {
            get => _Samples[n];
            set => _Samples[n] = value;
        }

        [NotNull] public double[] Samples => _Samples;

        public SamplesDigitalSignal(double dt, [NotNull] double[] Samples) : base(dt) => _Samples = Samples ?? throw new ArgumentNullException(nameof(Samples));

        public SamplesDigitalSignal(double dt, int SamplesCount, [NotNull] Func<double, double> f) : this(dt, f.Sampling(0, dt, SamplesCount)) { }

        public SamplesDigitalSignal(double dt, [NotNull] IEnumerable<double> Samples) : this(dt, (Samples ?? throw new ArgumentNullException(nameof(Samples))).ToArray()) { }
        public SamplesDigitalSignal(double dt, [NotNull] IEnumerable<int> Samples) : this(dt, (Samples ?? throw new ArgumentNullException(nameof(Samples))).Select(v => (double)v).ToArray()) { }

        [NotNull]
        public SamplesDigitalSignal GetIntegral(double s0 = 0)
        {
            var samples = new double[_Samples.Length];
            var dt05 = _dt / 2;
            samples[0] = s0;
            for (int i = 1, count = samples.Length; i < count; i++)
                samples[i] = samples[i - 1] + (_Samples[i - 1] + _Samples[i]) / dt05;

            return new SamplesDigitalSignal(_dt, samples);
        }

        #region Override Object

        public override string ToString() => $"signal dt:{_dt}; count:{_Samples.Length}; power:{Power.RoundAdaptive(2)}";

        public bool Equals(SamplesDigitalSignal other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            if (!_dt.Equals(other._dt) || _Samples.Length != other._Samples.Length) return false;
            if (ReferenceEquals(_Samples, other._Samples)) return true;
            for (var i = 0; i < _Samples.Length; i++)
                if (!_Samples[i].Equals(other._Samples[i]))
                    return false;
            return true;
        }

        public override bool Equals(object obj) =>
            obj != null
            && (ReferenceEquals(this, obj)
                || obj.GetType() == GetType() && Equals((SamplesDigitalSignal)obj));

        public override int GetHashCode()
        {
            unchecked
            {
                var hash = _dt.GetHashCode() * 397;
                for (var i = 0; i < _Samples.Length; i++)
                    hash = (hash * 397) ^ _Samples[i].GetHashCode();
                return hash;
            }
        }

        #endregion

        public override IEnumerator<double> GetEnumerator() => ((IEnumerable<double>) _Samples).GetEnumerator();

        [NotNull] public static implicit operator double[] ([NotNull] SamplesDigitalSignal s) => s._Samples;
    }
}
