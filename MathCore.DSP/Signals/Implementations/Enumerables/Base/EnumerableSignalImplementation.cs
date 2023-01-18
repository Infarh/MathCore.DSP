using System.Collections;

namespace MathCore.DSP.Signals.Implementations.Enumerables.Base;

public abstract class EnumerableSignalImplementation : EnumerableSignal
{
    protected EnumerableSignalImplementation(double dt, IEnumerable<double> Samples, double t0) : base(dt, Samples, t0) { }

    protected abstract class SignalInfo : IEnumerable<double>
    {
        protected double _A;
        protected double _dt;
        protected int _SamplesCount;

        public double A { get => _A; set => _A = value; }

        public double dt { get => _dt; set => _dt = value; }

        public int SamplesCount { get => _SamplesCount; set => _SamplesCount = value; }

        protected SignalInfo(double A, double dt, int SamplesCount) => (_A, _dt, _SamplesCount) = (A, dt, SamplesCount);

        protected abstract double Sample(double t);

        public virtual IEnumerator<double> GetEnumerator()
        {
            var count = _SamplesCount;
            //var infinity = count < 0;
            var t = 0d;
            var dt0 = _dt;
            while (count-- > 0) yield return Sample(t += dt0);
        }

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
    }
}