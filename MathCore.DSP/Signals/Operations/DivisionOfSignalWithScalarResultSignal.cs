using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    public class DivisionOfSignalWithScalarResultSignal : BinaryScalarOperationResultSignal
    {
        public DivisionOfSignalWithScalarResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => x / y) { }
    }
}