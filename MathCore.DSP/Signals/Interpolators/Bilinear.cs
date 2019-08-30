using System;
using System.Collections.Generic;
using System.Text;

namespace MathCore.DSP.Signals.Interpolators
{
    [Serializable]
    public class Bilinear : Interpolator
    {
        public override double GetValue(double[] Values, double pos) => throw new NotImplementedException();

        public override double SetValue(double[] Values, double pos, double Value) => throw new NotSupportedException();
    }
}
