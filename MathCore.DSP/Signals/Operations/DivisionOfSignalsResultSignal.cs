namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции деления сигнала на сигнал</summary>
public class DivisionOfSignalsResultSignal(DigitalSignal S1, DigitalSignal S2) : BinaryOperationResultSignal(S1.NotNull(), S2.NotNull(), (x, y) => x / y);