using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    /// <summary>Методы-расширения цифровых фильтров</summary>
    public static class DigitalFilterExtension
    {
        /// <summary>Получить импульсную характеристику фильтра</summary>
        /// <param name="filter">Объект фильтра, для которого требуется получить импульсную характеристику</param>
        /// <param name="MaxSamples">Максимальное количество отсчётов (если меньше 0, то число отсчётов ограничивается по точности)</param>
        /// <param name="Accuracity">Точность вычисления по мощности (энергии)</param>
        /// <returns>Последовательность отсчётов импульсной характеристики</returns>
        public static IEnumerable<double> GetImpulseResponse([NotNull] this DigitalFilter filter, int MaxSamples = -1, double Accuracity = 0.001)
        {
            if (filter is null) throw new ArgumentNullException(nameof(filter));
            if (Accuracity <= 0) throw new ArgumentOutOfRangeException(nameof(Accuracity), "Точность должна быть больше 0");
            if (MaxSamples == 0) yield break;

            var filter_order = filter.Order;
            var state = new double[filter_order];
            yield return filter.Process(1, state);
            var energy = state.Sum(s => s * s) / filter_order;
            var energy_max = energy;
            do
            {
                if (--MaxSamples == 0) yield break;
                yield return filter.Process(0, state);
                energy = state.Sum(s => s * s) / filter_order;
                if (energy > energy_max)
                    energy_max = energy;
            } while (energy / energy_max >= Accuracity);
        }
    }
}