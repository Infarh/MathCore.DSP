namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал, как результат операции умножения сигнала на число</summary>
public class MultiplyOfSignalWithScalarResultSignal : BinaryScalarOperationResultSignal
{
    public MultiplyOfSignalWithScalarResultSignal(DigitalSignal S, double X) : base(S.NotNull(), X, (x, y) => x * y) { }
}