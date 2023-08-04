namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции умножения сигнала на число</summary>
public class MultiplyOfSignalWithScalarResultSignal(DigitalSignal S, double X) : BinaryScalarOperationResultSignal(S.NotNull(), X, (x, y) => x * y);