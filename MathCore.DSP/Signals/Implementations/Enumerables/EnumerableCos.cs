using MathCore.DSP.Signals.Implementations.Enumerables.Base;

using static System.Math;

namespace MathCore.DSP.Signals.Implementations.Enumerables;

public class EnumerableCos(double A, double f0, double phi0, double dt, double t0 = 0, int SamplesCount = -1)
    : PeriodicSignal(dt, new CosSignalInfo(A, f0, phi0, dt, SamplesCount), t0)
{
    protected class CosSignalInfo(double A, double f0, double phi0, double dt, int SamplesCount)
        : PeriodicSignalInfo(A, f0, phi0, dt, SamplesCount)
    {
        protected override double Sample(double t) => _A * Cos(Consts.pi2 * _f0 * t + _phi0);
    }
}