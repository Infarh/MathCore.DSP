namespace MathCore.DSP.Signals.Operations;

public class SubstractionOfSignalWithScalarResultSignal(DigitalSignal S, double X) : BinaryScalarOperationResultSignal(S.NotNull(), X, (x, y) => x - y);