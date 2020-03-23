﻿using System;
using MathCore.DSP.Signals.Implementations.Enumerables.Base;

namespace MathCore.DSP.Signals.Implementations
{
    public class EnumerableCos : HarmonicSignal
    {
        protected class CosSignalCore : HarmonicSignalCore
        {
            public CosSignalCore(double A, double f0, double phi0, double dt, int SamplesCount) : base(A, f0, phi0, dt, SamplesCount) { }
            protected override double Sample(double t) => _A * Math.Cos(Math.PI * 2 * _f0 + _phi0);
        }

        public EnumerableCos(double A, double f0, double phi0, double dt, double t0 = 0, int SamplesCount = -1)
            : base(dt, new CosSignalCore(A, f0, phi0, dt, SamplesCount), t0) { }
    }
}