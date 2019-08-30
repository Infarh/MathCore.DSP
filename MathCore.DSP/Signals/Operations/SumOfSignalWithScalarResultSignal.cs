using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    public class SumOfSignalWithScalarResultSignal : BinaryScalarOperationResultSignal
    {
        public SumOfSignalWithScalarResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => x + y) { }
    }
}