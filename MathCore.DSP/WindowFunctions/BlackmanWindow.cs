using static System.Math;

namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно Блэкмана (Blackman)</summary>
/// <remarks>
/// Классическое косинусное окно с хорошим компромиссом между шириной главного лепестка и уровнем боковых лепестков; часто применяется при спектральном анализе и в проектировании FIR‑фильтров
/// Используется, когда важно снизить утечку спектра по сравнению с Хэннингом/Хэммингом при умеренной потере разрешения
/// Источники: R. B. Blackman, J. W. Tukey, The Measurement of Power Spectra, 1958; F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public class BlackmanWindow
{
    private const double alpha = 0.16; // классический параметр альфа

    private const double a0 = (1 - alpha) / 2; // коэффициент a0
    private const double a1 = alpha / 2;       // коэффициент a1

    private const double pi4 = Consts.pi2 * 2;

    /// <summary>Значение параметризованного окна Блэкмана в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <param name="a">Параметр альфа</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N, double a) => Value((double)n / N, a);

    /// <summary>Значение параметризованного окна Блэкмана в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <param name="a">Параметр альфа</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x, double a) => (1 - a) - Cos(Consts.pi2 * x) + a * Cos(pi4 * x);

    /// <summary>Значение стандартного окна Блэкмана (alpha = 0.16) в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N) => Value((double)n / N);

    /// <summary>Значение стандартного окна Блэкмана (alpha = 0.16) в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x) => a0 - 0.5 * Cos(Consts.pi2 * x) + a1 * Cos(pi4 * x);

    //public static double Value(double x) =>
    //    Abs(x) <= 0.5
    //        ? (25 * Cos(2 * PI * x) + 4 * Cos(4 * PI * x) + 21) / 50
    //        : 0;

    //public static double Nuttall(double x) =>
    //    Abs(x) <= 0.5
    //        ? (121_849 * Cos(2 * PI * x) + 36_058 * Cos(4 * PI * x) + 3_151 * Cos(6 * PI * x) + 88_942) / 250_000
    //        : 0;
}