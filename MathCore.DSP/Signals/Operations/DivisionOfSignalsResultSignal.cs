using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    /// <summary>Сигнал, как результат операции деления сигнала на сигнал</summary>
    public class DivisionOfSignalsResultSignal : BinaryOperationResultSignal
    {
        public DivisionOfSignalsResultSignal([NotNull] DigitalSignal S1, [NotNull] DigitalSignal S2) : base(S1, S2, (x, y) => x / y) { }
    }
}