namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно Кайзера (Kaiser)</summary>
/// <remarks>
/// Параметрическое окно на основе модифицированной функции Бесселя нулевого порядка; параметр alpha управляет компромиссом между шириной главного лепестка и уровнем боковых лепестков
/// Очень удобно при проектировании FIR‑фильтров (метод Кайзера) и при спектральном анализе, где нужен гибкий контроль характеристик
/// См.: J. F. Kaiser, Nonrecursive digital filter design using the I0-sinh window function, 1966; также F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class KaiserWindow
{
    /// <summary>Значение окна Кайзера в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <param name="alpha">Параметр окна (чем больше, тем уже главная полоса и ниже боковые лепестки)</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N, double alpha) => Value((double)n / N, alpha);

    /// <summary>Значение окна Кайзера в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <param name="alpha">Параметр окна (чем больше, тем уже главная полоса и ниже боковые лепестки)</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x, double alpha) => 
        SpecialFunctions.Bessel.I0(alpha * (1 - (2 * (x - 0.5)).Pow2()).Sqrt()) /
        SpecialFunctions.Bessel.I0(alpha);
}
