using MathCore.DSP.Signals.Implementations.Enumerables.Base;

using static System.Math;

namespace MathCore.DSP.Signals.Implementations.Enumerables;

public class EnumerableSin : PeriodicSignal
{
    protected class SinSignalInfo : PeriodicSignalInfo
    {
        public SinSignalInfo(double A, double f0, double phi0, double dt, int SamplesCount) : base(A, f0, phi0, dt, SamplesCount) { }
        protected override double Sample(double t) => _A * Sin(PI * 2 * _f0 * t + _phi0);
    }

    public EnumerableSin(double A, double f0, double phi0, double dt, double t0 = 0, int SamplesCount = -1) 
        : base(dt, new SinSignalInfo(A, f0, phi0, dt, SamplesCount), t0) { }
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