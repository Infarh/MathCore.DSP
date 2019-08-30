using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    public class SubstractionOfScalarWithSignalResultSignal : BinaryScalarOperationResultSignal
    {
        public SubstractionOfScalarWithSignalResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => y - x) { }
    }
}