using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    public class MultiplyOfSignalWithScalarResultSignal : BinaryScalarOperationResultSignal
    {
        public MultiplyOfSignalWithScalarResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => x * y) { }
    }
}