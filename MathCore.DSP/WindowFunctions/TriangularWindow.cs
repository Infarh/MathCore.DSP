namespace MathCore.DSP.WindowFunctions;

/// <summary>Треугольное окно (Triangular, Bartlett)</summary>
/// <remarks>
/// Линейно убывающее к краям окно с умеренной шириной главного лепестка и средним подавлением боковых лепестков
/// Подходит для простого спектрального анализа, фильтрации и сглаживания, когда требуется более мягкая аподизация, чем у прямоугольного окна
/// См.: F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class TriangularWindow
{
    /// <summary>Значение треугольного окна в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N) => Value((double)n / N);

    /// <summary>Значение треугольного окна в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x) => 1 - 2 * Math.Abs(x - 0.5);

    //public static double Value(int n, int N) => 1 - 2 * Math.Abs((double)n / N - 0.5) ;
    //public static double Value(int n, int N) => (N - 2 * Math.Abs(n - N / 2d)) / N;
}
