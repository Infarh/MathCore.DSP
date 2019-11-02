using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    /// <summary>Сигнал, как результат операции вычитания сигнала из сигнала</summary>
    public class SubstractionOfSignalsResultSignal : BinaryOperationResultSignal
    {
        public SubstractionOfSignalsResultSignal([NotNull] DigitalSignal S1, [NotNull] DigitalSignal S2) : base(S1, S2, (x, y) => x - y) { }
    }
}