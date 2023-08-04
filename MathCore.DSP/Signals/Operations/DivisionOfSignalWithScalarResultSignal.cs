namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат бинарной операции между сигналом и числом</summary>
public class DivisionOfSignalWithScalarResultSignal(DigitalSignal S, double X) : BinaryScalarOperationResultSignal(S.NotNull(), X, (x, y) => x / y);