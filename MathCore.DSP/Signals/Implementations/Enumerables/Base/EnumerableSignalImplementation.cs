using System.Collections;

namespace MathCore.DSP.Signals.Implementations.Enumerables.Base;

public abstract class EnumerableSignalImplementation(double dt, IEnumerable<double> Samples, double t0) : EnumerableSignal(dt, Samples, t0)
{
    protected abstract class SignalInfo(double A, double dt, int SamplesCount) : IEnumerable<double>
    {
        protected double _A = A;
        protected double _dt = dt;
        protected int _SamplesCount = SamplesCount;

        public double A { get => _A; set => _A = value; }

        public double dt { get => _dt; set => _dt = value; }

        public int SamplesCount { get => _SamplesCount; set => _SamplesCount = value; }

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