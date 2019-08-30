using System;

namespace MathCore.DSP.Filters
{
    public class BandPassRLC : IIR
    {
        public BandPassRLC(double f0, double DeltaF, double dt)
            : this(Math.Tan(Math.PI * f0 * dt), Math.PI * DeltaF * dt)
        { }

        private BandPassRLC(double w0, double dw)
            : base(
                A: new[] { w0 * w0 + dw + 1, 2 * (w0 * w0 - 1), w0 * w0 - dw + 1 },
                B: new[] { dw, 0, -dw })
        { }
    }
}