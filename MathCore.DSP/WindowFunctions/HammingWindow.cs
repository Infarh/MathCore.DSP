namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно Хэмминга (Hamming)</summary>
/// <remarks>
/// Косинусное окно, снижающее утечку спектра за счёт уменьшения боковых лепестков при небольшой ширине главного лепестка; применимо в спектральном анализе и при проектировании FIR‑фильтров
/// Выбирается как компромисс между окном Ханна и прямоугольным, когда нужно умеренное подавление утечки и неплохое частотное разрешение
/// Источники: R. W. Hamming, Error Detecting and Error Correcting Codes, 1950; F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public class HammingWindow
{
    /// <summary>Значение параметризованного окна Хэмминга в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <param name="a0">Коэффициент a0 (обычно ≈ 0.54)</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N, double a0) => Value((double)n / N, a0);

    /// <summary>Значение параметризованного окна Хэмминга в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <param name="a0">Коэффициент a0 (обычно ≈ 0.54)</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x, double a0) => a0 - (1 - a0) * Math.Cos(Consts.pi2 * x);

    private const double a0 = 25d / 46; // стандартный коэффициент a0 ≈ 0.543478
    private const double a1 = 1 - a0;   // вспомогательный коэффициент

    /// <summary>Значение стандартного окна Хэмминга в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N) => Value((double)n / N);

    /// <summary>Значение стандартного окна Хэмминга в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x) => a0 - a1 * Math.Cos(Consts.pi2 * x);
}
