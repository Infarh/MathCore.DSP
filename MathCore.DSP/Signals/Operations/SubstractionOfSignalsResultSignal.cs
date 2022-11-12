namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции вычитания сигнала из сигнала</summary>
public class SubstractionOfSignalsResultSignal : BinaryOperationResultSignal
{
    public SubstractionOfSignalsResultSignal(DigitalSignal S1, DigitalSignal S2) : base(S1.NotNull(), S2.NotNull(), (x, y) => x - y) { }
}