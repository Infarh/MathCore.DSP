// ReSharper disable InconsistentNaming
namespace MathCore.DSP.Signals.Implementations.Enumerables.Base
{
    public abstract class HarmonicSignal : EnumerableSignalImplementation
    {
        protected abstract class HarmonicSignalCore : SignalCore
        {
            protected double _f0;
            protected double _phi0;

            public double f0 { get => _f0; set => _f0 = value; }
            public double phi0 { get => _phi0; set => _phi0 = value; }

            protected HarmonicSignalCore(double A, double f0, double phi0, double dt, int SamplesCount) : base(A, dt, SamplesCount)
            {
                _f0 = f0;
                _phi0 = phi0;
            }
        }

        private readonly HarmonicSignalCore _SignalCore;

        public double A { get => _SignalCore.A; set => _SignalCore.A = value; }

        public double f0 { get => _SignalCore.f0; set => _SignalCore.f0 = value; }

        public double phi0 { get => _SignalCore.phi0; set => _SignalCore.phi0 = value; }
            
        protected HarmonicSignal(double dt, HarmonicSignalCore SignalCore, double t0) : base(dt, SignalCore, t0) => _SignalCore = SignalCore;
    }
}