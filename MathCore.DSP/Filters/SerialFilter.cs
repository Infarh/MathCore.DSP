using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    /// <summary>Последовательный комбинационный фильтр</summary>
    public class SerialFilter : CombinationFilter
    {
        /// <summary>Инициализация нового последовательного комбинационного фильтра</summary>
        /// <param name="Filter1">Первый фильтр в комбинации</param>
        /// <param name="Filter2">Второй фильтр в комбинации</param>
        public SerialFilter([NotNull] Filter Filter1, [NotNull] Filter Filter2) : base(Filter1, Filter2) { }

        public override double Process(double Sample) => Filter2.Process(Filter1.Process(Sample));

        public override void Reset()
        {
            Filter1.Reset();
            Filter2.Reset();
        }

        public override Complex GetTransmissionCoefficient(double f) => Filter1.GetTransmissionCoefficient(f) * Filter2.GetTransmissionCoefficient(f);
    }
}