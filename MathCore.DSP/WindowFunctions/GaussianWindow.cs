namespace MathCore.DSP.WindowFunctions;

/// <summary>Гауссово окно (Gaussian)</summary>
/// <remarks>
/// Гладкое окно на основе нормального распределения, обеспечивает хорошее подавление дальних боковых лепестков при управляемой ширине через параметр alpha
/// Подходит для задач, где требуется минимизация ряби и хорошее временно-частотное локализующее свойство (например, в обработке аудиосигналов)
/// См.: F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class GaussianWindow
{
    /// <summary>Значение Гауссова окна в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <param name="alpha">Параметр ширины окна (чем больше, тем шире)</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N, double alpha) => Value((double)n / N, alpha);

    /// <summary>Значение Гауссова окна в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <param name="alpha">Параметр ширины окна (чем больше, тем шире)</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x, double alpha) => Math.Exp(-2 * ((x - 0.5) / alpha).Pow2());
}
