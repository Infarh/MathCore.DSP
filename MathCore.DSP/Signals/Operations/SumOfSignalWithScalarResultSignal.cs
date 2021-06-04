namespace MathCore.DSP.Signals.Operations
{
    /// <summary>Сигнал, как результат операции сложения сигнала с числом</summary>
    public class SumOfSignalWithScalarResultSignal : BinaryScalarOperationResultSignal
    {
        public SumOfSignalWithScalarResultSignal(DigitalSignal S, double X) : base(S, X, (x, y) => x + y) { }
    }
}