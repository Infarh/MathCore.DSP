namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно с плоской вершиной (Flat Top)</summary>
/// <remarks>
/// Косинусное окно, обеспечивающее очень плоскую вершину амплитудной характеристики, что уменьшает погрешность амплитудных измерений при спектральном анализе
/// Рекомендуется для точного измерения амплитуд гармоник (например, калибровка и измерительная FFT), когда небольшое расширение главного лепестка допустимо
/// См.: IEEE Std 1057 и 1241; F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class FlatTopWindow
{
    private const double a1 = 1.93;
    private const double a2 = 1.29;
    private const double a3 = 0.388;
    private const double a4 = 0.032;

    private const double pi2 = Consts.pi2 * 1;
    private const double pi4 = Consts.pi2 * 2;
    private const double pi6 = Consts.pi2 * 3;
    private const double pi8 = Consts.pi2 * 4;

    /// <summary>Значение окна в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N) => Value((double)n / N);

    /// <summary>Значение окна в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x) =>
        1
        - a1 * Math.Cos(pi2 * x)
        + a2 * Math.Cos(pi4 * x)
        - a3 * Math.Cos(pi6 * x)
        + a4 * Math.Cos(pi8 * x);

    /// <summary>Значение параметризованного окна Flat Top в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <param name="a1">Коэффициент a1</param>
    /// <param name="a2">Коэффициент a2</param>
    /// <param name="a3">Коэффициент a3</param>
    /// <param name="a4">Коэффициент a4</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N, double a1, double a2, double a3, double a4) => Value((double)n / N, a1, a2, a3, a4);

    /// <summary>Значение параметризованного окна Flat Top в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <param name="a1">Коэффициент a1</param>
    /// <param name="a2">Коэффициент a2</param>
    /// <param name="a3">Коэффициент a3</param>
    /// <param name="a4">Коэффициент a4</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x, double a1, double a2, double a3, double a4) =>
        1
        - a1 * Math.Cos(pi2 * x)
        + a2 * Math.Cos(pi4 * x)
        - a3 * Math.Cos(pi6 * x)
        + a4 * Math.Cos(pi8 * x);
}
