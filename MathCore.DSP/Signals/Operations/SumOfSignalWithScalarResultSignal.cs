namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции сложения сигнала с числом</summary>
public class SumOfSignalWithScalarResultSignal(DigitalSignal S, double X) : BinaryScalarOperationResultSignal(S.NotNull(), X, (x, y) => x + y);