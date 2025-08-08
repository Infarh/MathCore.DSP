namespace MathCore.DSP.Samples.Extensions;

public static class SampleSI16Ex
{
    /// <summary>Фазовая демодуляция радиосигнала</summary>
    /// <param name="samples">Последовательность отсчётов квадратурного радиосигнала</param>
    /// <param name="f0">Центральная частота фазовой модуляции</param>
    /// <param name="fd">Частота дискретизации</param>
    /// <returns>Возвращает массив вещественных значений отсчётов демодулированного сигнала</returns>
    public static float[] PhaseDemodulation(this Span<SampleSI16> samples, double f0, double fd)
    {
        throw new NotImplementedException();
    }
}
