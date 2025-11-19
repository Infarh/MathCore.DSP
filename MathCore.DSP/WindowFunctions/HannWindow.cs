namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно Ханна (Hann, Хэннинг)</summary>
/// <remarks>
/// Простое косинусное окно с хорошим компромиссом между разрешением и подавлением утечки; одно из самых распространённых для общего спектрального анализа
/// Рекомендуется как базовый выбор для FFT‑анализаторов, когда требуется умеренная ширина главного лепестка и пониженные боковые лепестки
/// См.: J. W. Tukey, 1967; F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class HannWindow
{
    /// <summary>Значение окна Ханна в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N) => Value((double)n / N);

    /// <summary>Значение окна Ханна в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x) => Math.Sin(Math.PI * x).Pow2();
    //public static double Value(int n, int N) => 0.5 - 0.5 * Math.Cos(Math.PI * n / N);
}
