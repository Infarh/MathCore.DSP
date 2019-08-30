using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    public static class DigitalFilterExtension
    {
        public static IEnumerable<double> GetImpulseResponse([NotNull] this DigitalFilter filter, int MaxSamples = -1, double accuracity = 0.001)
        {
            if (filter is null) throw new ArgumentNullException(nameof(filter));
            if (accuracity <= 0) throw new ArgumentOutOfRangeException(nameof(accuracity), "Точность должна быть больше 0");
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
            } while (energy / energy_max >= accuracity);
        }
    }
}