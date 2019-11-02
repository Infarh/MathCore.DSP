using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    /// <summary>Сигнал, как результат операции умножения сигнала на число</summary>
    public class MultiplyOfSignalWithScalarResultSignal : BinaryScalarOperationResultSignal
    {
        public MultiplyOfSignalWithScalarResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => x * y) { }
    }
}