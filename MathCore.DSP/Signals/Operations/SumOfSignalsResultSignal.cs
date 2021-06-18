namespace MathCore.DSP.Signals.Operations
{
    /// <summary>Сигнал, как результат операции сложения сигналов</summary>
    public class SumOfSignalsResultSignal : BinaryOperationResultSignal
    {
        public SumOfSignalsResultSignal(DigitalSignal S1, DigitalSignal S2) : base(S1, S2, (x, y) => x + y) { }
    }
}