using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    /// <summary>Сигнал, как результат операции деления сигнала на число</summary>
    public class DivisionOfScalarWithSignalResultSignal : BinaryScalarOperationResultSignal
    {
        public DivisionOfScalarWithSignalResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => x / y) { }
    }
}