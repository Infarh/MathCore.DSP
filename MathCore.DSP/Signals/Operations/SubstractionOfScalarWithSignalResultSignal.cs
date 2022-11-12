namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции вычитания числа из сигнала</summary>
public class SubstractionOfScalarWithSignalResultSignal : BinaryScalarOperationResultSignal
{
    public SubstractionOfScalarWithSignalResultSignal(DigitalSignal S, double X) : base(S.NotNull(), X, (x, y) => y - x) { }
}