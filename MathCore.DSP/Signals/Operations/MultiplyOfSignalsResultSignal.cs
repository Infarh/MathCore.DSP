namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции умножения двух сигналов</summary>
public class MultiplyOfSignalsResultSignal(DigitalSignal S1, DigitalSignal S2) : BinaryOperationResultSignal(S1.NotNull(), S2.NotNull(), (x, y) => x * y);