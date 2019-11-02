using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    /// <summary>Сигнал, как результат бинарной операции между сигналом и числом</summary>
    public class DivisionOfSignalWithScalarResultSignal : BinaryScalarOperationResultSignal
    {
        public DivisionOfSignalWithScalarResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => x / y) { }
    }
}