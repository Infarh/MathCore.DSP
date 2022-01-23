namespace MathCore.DSP.Signals.Interpolators;

/// <summary>Интерполятор</summary>
[Serializable]
public class Interpolator : ICloneable
{
    public double this[double pos, double[] Values] { get => GetValue(Values, pos); set => SetValue(Values, pos, value); }

    /// <summary>Получить значение</summary>
    /// <param name="Values">Массив значений</param>
    /// <param name="pos">Положение в массиве</param>
    /// <returns>Интерполированное значение</returns>
    public virtual double GetValue(double[] Values, double pos) => Values[(int)pos];

    /// <summary>Установить значение</summary>
    /// <param name="Values">Массив значений</param>
    /// <param name="pos">Положение в массиве</param>
    /// <param name="Value">Устанавливаемое значение</param>
    /// <returns>Установленное значение</returns>
    public virtual double SetValue(double[] Values, double pos, double Value) => Values[(int)pos] = Value;

    #region ICloneable Members

    public virtual object Clone() => MemberwiseClone();

    #endregion
}