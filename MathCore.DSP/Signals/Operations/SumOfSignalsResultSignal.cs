using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    public class SumOfSignalsResultSignal : BinaryOperationResultSignal
    {
        public SumOfSignalsResultSignal([NotNull] DigitalSignal S1, [NotNull] DigitalSignal S2) : base(S1, S2, (x, y) => x + y) { }
    }
}