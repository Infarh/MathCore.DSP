namespace MathCore.DSP.Signals.Operations
{
    /// <summary>Сигнал - результат операции над сигналами</summary>
    public abstract class OperationResultSignal : DigitalSignal
    {
        /// <summary>Инициализация нового результирующего сигнала</summary>
        /// <param name="dt">Период дискретизации</param>
        protected OperationResultSignal(double dt) : base(dt) { }
    }
}