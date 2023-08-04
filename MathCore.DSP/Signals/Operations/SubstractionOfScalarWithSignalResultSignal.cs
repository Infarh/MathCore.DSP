namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции вычитания числа из сигнала</summary>
public class SubstractionOfScalarWithSignalResultSignal(DigitalSignal S, double X) : BinaryScalarOperationResultSignal(S.NotNull(), X, (x, y) => y - x);