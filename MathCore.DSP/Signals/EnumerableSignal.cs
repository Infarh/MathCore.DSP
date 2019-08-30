using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;

namespace MathCore.DSP.Signals
{
    public class EnumerableSignal : DigitalSignal
    {
        private readonly IEnumerable<double> _Samples;
        public override int SamplesCount => _Samples.Count();

        public override double this[int n]
        {
            get => throw new NotSupportedException();
            set => throw new NotSupportedException();
        }

        public EnumerableSignal(double dt, [NotNull] IEnumerable<double> Samples) : base(dt) => _Samples = Samples ?? throw new ArgumentNullException(nameof(Samples));

        public override IEnumerator<double> GetEnumerator() => _Samples.GetEnumerator();
    }
}