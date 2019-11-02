using System;
using MathCore.Annotations;

namespace MathCore.DSP.Signals.Interpolators
{
    /// <summary>Интерполятор</summary>
    [Serializable]
    public class Interpolator : ICloneable
    {
        public double this[double pos, [NotNull] double[] Values] { get => GetValue(Values, pos); set => SetValue(Values, pos, value); }

        /// <summary>Получить значение</summary>
        /// <param name="Values">Массив значений</param>
        /// <param name="pos">Положение в массиве</param>
        /// <returns>Интерполированное значение</returns>
        public virtual double GetValue([NotNull] double[] Values, double pos) => Values[(int)pos];

        /// <summary>Установить значение</summary>
        /// <param name="Values">Массив значений</param>
        /// <param name="pos">Положение в массиве</param>
        /// <param name="Value">Устанавливаемое значение</param>
        /// <returns>Установленное значение</returns>
        public virtual double SetValue([NotNull] double[] Values, double pos, double Value) => Values[(int)pos] = Value;

        #region ICloneable Members

        public virtual object Clone() => MemberwiseClone();

        #endregion
    }
}
