using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;
using MathCore.DSP.Signals;

namespace MathCore.DSP.Filters
{
    public abstract class Filter
    {
        public abstract double Process(double sample);

        [NotNull] public virtual IEnumerable<double> Process([NotNull] IEnumerable<double> samples) => samples.Select(Process);

        [NotNull] public DigitalSignal Process([NotNull] DigitalSignal s) => new SamplesDigitalSignal((s ?? throw new ArgumentNullException(nameof(s))).dt, Process((IEnumerable<double>)s));

        public abstract void Reset();

        public abstract Complex GetTransmissionCoefficient(double f);

        [NotNull]
        public static DigitalSignal operator *([NotNull] Filter filter, [NotNull] DigitalSignal signal)
        {
            if (filter is null) throw new ArgumentNullException(nameof(filter));
            if (signal is null) throw new ArgumentNullException(nameof(signal));
            return filter.Process(signal);
        }

        [NotNull] public static DigitalSignal operator *([NotNull] DigitalSignal signal, [NotNull] Filter filter) => filter * signal;

        [NotNull] public static ParallelFilter operator +([NotNull] Filter filter1, [NotNull] Filter filter2) => new ParallelFilter(filter1, filter2);

        [NotNull] public static SerialFilter operator *([NotNull] Filter filter1, [NotNull] Filter filter2) => new SerialFilter(filter1, filter2);
    }

    public abstract class CombinationFilter : Filter
    {
        [NotNull] public Filter Filter1 { get; }
        [NotNull] public Filter Filter2 { get; }

        protected CombinationFilter([NotNull] Filter Filter1, [NotNull] Filter Filter2)
        {
            this.Filter1 = Filter1 ?? throw new ArgumentNullException(nameof(Filter1));
            this.Filter2 = Filter2 ?? throw new ArgumentNullException(nameof(Filter2));
        }
    }

    public class ParallelFilter : CombinationFilter
    {

        public ParallelFilter(Filter Filter1, Filter Filter2) : base(Filter1, Filter2) { }

        public override double Process(double sample) => Filter1.Process(sample / 2) + Filter2.Process(sample / 2);

        public override void Reset()
        {
            Filter1.Reset();
            Filter2.Reset();
        }

        public override Complex GetTransmissionCoefficient(double f) => (Filter1.GetTransmissionCoefficient(f) + Filter2.GetTransmissionCoefficient(f)) / 2;
    }

    public class SerialFilter : CombinationFilter
    {

        public SerialFilter(Filter Filter1, Filter Filter2) : base(Filter1, Filter2) { }

        public override double Process(double sample) => Filter2.Process(Filter1.Process(sample));

        public override void Reset()
        {
            Filter1.Reset();
            Filter2.Reset();
        }

        public override Complex GetTransmissionCoefficient(double f) => Filter1.GetTransmissionCoefficient(f) * Filter2.GetTransmissionCoefficient(f);
    }
}