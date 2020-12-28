using System;

namespace MathCore.DSP.Signals.Interpolators
{
    /// <summary>Билинейный интерполятор</summary> 
    [Serializable]
    public class Bilinear : Interpolator
    {
        public override double GetValue(double[] Values, double pos) => throw new NotImplementedException();

        public override double SetValue(double[] Values, double pos, double Value) => throw new NotSupportedException();
    }
}
