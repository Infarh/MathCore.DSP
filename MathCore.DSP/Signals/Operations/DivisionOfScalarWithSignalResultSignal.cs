namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции деления сигнала на число</summary>
public class DivisionOfScalarWithSignalResultSignal(DigitalSignal S, double X) : BinaryScalarOperationResultSignal(S.NotNull(), X, (x, y) => x / y);