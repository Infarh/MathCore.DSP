namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно Бартлетта–Ханна (Bartlett–Hann)</summary>
/// <remarks>
/// Это окно сочетает свойства треугольного и косинусного окон и обеспечивает компромисс между шириной главного лепестка и уровнем боковых лепестков
/// Рекомендуется для спектрального анализа, когда требуется умеренная разрешающая способность и сниженная утечка спектра по сравнению с прямоугольным/треугольным окнами
/// Подробный обзор оконных функций и их свойств приведён в: F. J. Harris, On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform, Proceedings of the IEEE, 1978, DOI: 10.1109/PROC.1978.10837, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class BartlettHannWindow
{
    /// <summary>Значение окна в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N) => Value((double)n / N);

    /// <summary>Значение окна в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x) => 0.62 - 0.48 * Math.Abs(x - 0.5) - 0.38 * Math.Cos(Consts.pi2 * x);
}
