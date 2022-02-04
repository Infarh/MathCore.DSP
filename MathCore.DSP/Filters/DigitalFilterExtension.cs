using MathCore.Annotations;

namespace MathCore.DSP.Filters;

/// <summary>Методы-расширения цифровых фильтров</summary>
public static class DigitalFilterExtension
{
    /// <summary>Получить импульсную характеристику фильтра</summary>
    /// <param name="filter">Объект фильтра, для которого требуется получить импульсную характеристику</param>
    /// <param name="MaxSamples">Максимальное количество отсчётов (если меньше 0, то число отсчётов ограничивается по точности)</param>
    /// <param name="Accuracy">Точность вычисления по мощности (энергии) состояния фильтра</param>
    /// <returns>Последовательность отсчётов импульсной характеристики</returns>
    /// <exception cref="ArgumentOutOfRangeException">Если указанная точность меньше, либо равна 0</exception>
    public static IEnumerable<double> GetImpulseResponse([NotNull] this DigitalFilter filter, int MaxSamples = -1, double Accuracy = 0.001)
    {
        if (filter is null) throw new ArgumentNullException(nameof(filter));
        if (Accuracy <= 0) throw new ArgumentOutOfRangeException(nameof(Accuracy), "Точность должна быть больше 0");
        if (MaxSamples == 0) yield break;

        var state_length = filter.Order + 1;
        var state = new double[state_length];
        yield return filter.Process(1, state);
        var energy = state.Sum(s => s * s) / state_length;
        var energy_max = energy;
        do
        {
            if (--MaxSamples == 0) yield break;
            yield return filter.Process(0, state);
            energy = state.Sum(s => s * s) / state_length;
            if (energy > energy_max)
                energy_max = energy;
        } while (energy / energy_max >= Accuracy);
    }

    // <param name="Accuracy">Точность вычисления по мощности (энергии) состояния фильтра</param>
    /// <summary>Получить переходную характеристику фильтра</summary>
    /// <param name="filter">Объект фильтра, для которого требуется получить переходную характеристику</param>
    /// <param name="MaxSamples">Максимальное количество отсчётов (если меньше 0, то число отсчётов ограничивается по точности)</param>
    /// <returns>Последовательность отсчётов импульсной характеристики</returns>
    public static IEnumerable<double> GetTransientResponse([NotNull] this DigitalFilter filter, int MaxSamples = -1/*, double Accuracy = 0.001*/)
    {
        if (filter is null) throw new ArgumentNullException(nameof(filter));
        //if (Accuracy <= 0) throw new ArgumentOutOfRangeException(nameof(Accuracy), "Точность должна быть больше 0");
        if (MaxSamples == 0) yield break;

        var state_length = filter.Order + 1;
        var state = new double[state_length];
        yield return filter.Process(1, state);
        //var energy = state.Sum(s => s * s) / state_length;
        //var energy_max = energy;
        do
        {
            if (--MaxSamples == 0) yield break;
            yield return filter.Process(0, state);
            //energy = state.Sum(s => s * s) / state_length;
            //if (energy > energy_max)
            //    energy_max = energy;
        } while (/*energy / energy_max >= Accuracy*/ true);
    }
}