namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции деления сигнала на сигнал</summary>
public class DivisionOfSignalsResultSignal : BinaryOperationResultSignal
{
    public DivisionOfSignalsResultSignal(DigitalSignal S1, DigitalSignal S2) : base(S1, S2, (x, y) => x / y) { }
}