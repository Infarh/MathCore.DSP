namespace MathCore.DSP.Signals.Operations;

public class SubstractionOfSignalWithScalarResultSignal : BinaryScalarOperationResultSignal
{
    public SubstractionOfSignalWithScalarResultSignal(DigitalSignal S, double X) : base(S, X, (x, y) => x - y) { }
}