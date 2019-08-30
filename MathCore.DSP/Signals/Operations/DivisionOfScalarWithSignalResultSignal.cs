using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    public class DivisionOfScalarWithSignalResultSignal : BinaryScalarOperationResultSignal
    {
        public DivisionOfScalarWithSignalResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => x / y) { }
    }
}