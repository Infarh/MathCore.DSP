using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    public class SubstractionOfSignalWithScalarResultSignal : BinaryScalarOperationResultSignal
    {
        public SubstractionOfSignalWithScalarResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => x - y) { }
    }
}