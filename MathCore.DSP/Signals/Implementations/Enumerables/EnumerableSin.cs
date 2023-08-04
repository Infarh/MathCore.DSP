using MathCore.DSP.Signals.Implementations.Enumerables.Base;

using static System.Math;

namespace MathCore.DSP.Signals.Implementations.Enumerables;

public class EnumerableSin(double A, double f0, double phi0, double dt, double t0 = 0, int SamplesCount = -1)
    : PeriodicSignal(dt, new SinSignalInfo(A, f0, phi0, dt, SamplesCount), t0)
{
    protected class SinSignalInfo(double A, double f0, double phi0, double dt, int SamplesCount)
        : PeriodicSignalInfo(A, f0, phi0, dt, SamplesCount)
    {
        protected override double Sample(double t) => _A * Sin(Consts.pi2 * _f0 * t + _phi0);
    }
}

//public class EnumerableRect : PeriodicSignal
//{
//    protected class RectSignalInfo : PeriodicSignalInfo
//    {
//        public RectSignalInfo(double A, double f0, double phi0, double dt, int SamplesCount) : base(A, f0, phi0, dt, SamplesCount) { }
//        protected override double Sample(double t) => (t % T0);
//    }

//    public EnumerableRect(double A, double f0, double phi0, double dt, double t0 = 0, int SamplesCount = -1)
//        : base(dt, new RectSignalInfo(A, f0, phi0, dt, SamplesCount), t0) { }
//}