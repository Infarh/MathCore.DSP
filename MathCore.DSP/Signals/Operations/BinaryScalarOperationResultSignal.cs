using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    public abstract class BinaryScalarOperationResultSignal : OperationResultSignal
    {
        private readonly double _Value;
        [NotNull] private readonly Func<double, double, double> _Function;
        [NotNull] public DigitalSignal Signal { get; }

        public override int SamplesCount => Signal.SamplesCount;

        public override double this[int n]
        {
            get => _Function(Signal[n], _Value);
            set => throw new NotSupportedException();
        }

        protected BinaryScalarOperationResultSignal([NotNull] DigitalSignal Signal, double Value, [NotNull] Func<double, double, double> Function) : base(Signal.dt)
        {
            this.Signal = Signal ?? throw new ArgumentNullException(nameof(Signal));
            _Value = Value;
            _Function = Function ?? throw new ArgumentNullException(nameof(Function));
        }

        public override IEnumerator<double> GetEnumerator() => Signal.Select(s => _Function(s, _Value)).GetEnumerator();
    }
}