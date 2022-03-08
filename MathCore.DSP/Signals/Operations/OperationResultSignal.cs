namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал - результат операции над сигналами</summary>
public abstract class OperationResultSignal : DigitalSignal
{
    /// <summary>Инициализация нового результирующего сигнала</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="t0">Смещение во времени</param>
    protected OperationResultSignal(double dt, double t0 = 0) : base(dt, t0) { }
}