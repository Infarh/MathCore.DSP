﻿using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    /// <summary>Сигнал, как результат операции вычитания числа из сигнала</summary>
    public class SubstractionOfScalarWithSignalResultSignal : BinaryScalarOperationResultSignal
    {
        public SubstractionOfScalarWithSignalResultSignal([NotNull] DigitalSignal S, double X) : base(S, X, (x, y) => y - x) { }
    }
}