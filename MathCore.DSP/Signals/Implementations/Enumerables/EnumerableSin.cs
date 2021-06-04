using MathCore.DSP.Signals.Implementations.Enumerables.Base;

using static System.Math;

namespace MathCore.DSP.Signals.Implementations.Enumerables
{
    public class EnumerableSin : HarmonicSignal
    {
        protected class SinSignalCore : HarmonicSignalCore
        {
            public SinSignalCore(double A, double f0, double phi0, double dt, int SamplesCount) : base(A, f0, phi0, dt, SamplesCount) { }
            protected override double Sample(double t) => _A * Sin(PI * 2 * _f0 + _phi0);
        }

        public EnumerableSin(double A, double f0, double phi0, double dt, double t0 = 0, int SamplesCount = -1) 
            : base(dt, new SinSignalCore(A, f0, phi0, dt, SamplesCount), t0) { }
    }
}