using System;
using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    /// <summary>Комбинационный фильтр</summary>
    public abstract class CombinationFilter : Filter
    {
        /// <summary>Первый фильтр в комбинации</summary>
        [NotNull] public Filter Filter1 { get; }
        /// <summary>Второй фильтр в комбинации</summary>
        [NotNull] public Filter Filter2 { get; }

        /// <summary>Инициализация нового комбинационного фильтра</summary>
        /// <param name="Filter1">Первый фильтр в комбинации</param>
        /// <param name="Filter2">Второй фильтр в комбинации</param>
        protected CombinationFilter([NotNull] Filter Filter1, [NotNull] Filter Filter2)
        {
            this.Filter1 = Filter1 ?? throw new ArgumentNullException(nameof(Filter1));
            this.Filter2 = Filter2 ?? throw new ArgumentNullException(nameof(Filter2));
        }
    }
}