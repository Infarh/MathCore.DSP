namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции деления сигнала на число</summary>
public class DivisionOfScalarWithSignalResultSignal : BinaryScalarOperationResultSignal
{
    public DivisionOfScalarWithSignalResultSignal(DigitalSignal S, double X) : base(S.NotNull(), X, (x, y) => x / y) { }
}