using MathCore.DSP.Signals.Implementations.Enumerables.Base;

using static System.Math;

namespace MathCore.DSP.Signals.Implementations.Enumerables
{
    public class EnumerableCos : PeriodicSignal
    {
        protected class CosSignalInfo : PeriodicSignalInfo
        {
            public CosSignalInfo(double A, double f0, double phi0, double dt, int SamplesCount) : base(A, f0, phi0, dt, SamplesCount) { }
            protected override double Sample(double t) => _A * Cos(PI * 2 * _f0 * t + _phi0);
        }

        public EnumerableCos(double A, double f0, double phi0, double dt, double t0 = 0, int SamplesCount = -1)
            : base(dt, new CosSignalInfo(A, f0, phi0, dt, SamplesCount), t0) { }
    }
}