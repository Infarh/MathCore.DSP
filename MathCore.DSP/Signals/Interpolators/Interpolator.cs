using System;

namespace MathCore.DSP.Signals.Interpolators
{
    [Serializable]
    public class Interpolator : ICloneable
    {
        public double this[double pos, double[] Values] { get => GetValue(Values, pos); set => SetValue(Values, pos, value); }

        public virtual double GetValue(double[] Values, double pos) => Values[(int)pos];

        public virtual double SetValue(double[] Values, double pos, double Value) => Values[(int)pos] = Value;

        #region ICloneable Members

        public virtual object Clone() => MemberwiseClone();

        #endregion
    }
}
