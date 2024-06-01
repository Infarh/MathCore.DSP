namespace MathCore.DSP.Signals;

/// <summary>Класс методов-расширений для отсчётов спектра</summary>
public static class SpectrumSampleEx
{
    /// <summary>Представление последовательности отсчётов в виде последовательности комплексных чисел - значений</summary>
    /// <param name="samples">Последовательность отсчётов спектра</param>
    /// <returns>Последовательность комплексных чисел - отсчётов спектра</returns>
    public static IEnumerable<Complex> AsComplex(this IEnumerable<SpectrumSample> samples) => samples.Select(s => s.Value);

    /// <summary>Представление последовательности комплексных чисел в виде последовательности отсчётов спектра</summary>
    /// <param name="samples">Последовательность комплексных чисел - значений отсчётов спектра</param>
    /// <param name="df">Разрешающая способность спектра</param>
    /// <param name="f0">Начальное значение частоты - смещение</param>
    /// <returns>Последовательность отсчётов спектра</returns>
    public static IEnumerable<SpectrumSample> AsSamples(this IEnumerable<Complex> samples, double df, double f0 = 0)
    {
        var index = 0;
        foreach (var sample in samples)
        {
            yield return new(index * df + f0, sample);
            index++;
        }
    }
}