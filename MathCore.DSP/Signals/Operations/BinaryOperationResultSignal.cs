using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;

namespace MathCore.DSP.Signals.Operations
{
    public abstract class BinaryOperationResultSignal : OperationResultSignal
    {
        [NotNull] private readonly Func<double, double, double> _Function;
        [NotNull] public DigitalSignal Signal1 { get; }
        [NotNull] public DigitalSignal Signal2 { get; }

        public override int SamplesCount => Signal1.SamplesCount;

        public override double this[int n]
        {
            get => _Function(Signal1[n], Signal2[n]);
            set => throw new NotSupportedException();
        }

        protected BinaryOperationResultSignal([NotNull] DigitalSignal Signal1, [NotNull] DigitalSignal Signal2, [NotNull] Func<double, double, double> Function) : base(Signal1.dt)
        {
            this.Signal1 = Signal1 ?? throw new ArgumentNullException(nameof(Signal1));
            this.Signal2 = Signal2 ?? throw new ArgumentNullException(nameof(Signal2));
            _Function = Function ?? throw new ArgumentNullException(nameof(Function));
        }

        public override IEnumerator<double> GetEnumerator() => Signal1.Zip(Signal2, _Function).GetEnumerator();
    }
}