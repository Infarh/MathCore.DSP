namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно Блэкмана–Харриса (Blackman–Harris)</summary>
/// <remarks>
/// Многочленное окно с очень низкими боковыми лепестками при умеренной ширине главного лепестка, часто используется для высокоточного спектрального анализа амплитуд
/// Подходит, когда приоритетно подавление утечки спектра и минимум боковых лепестков ценой некоторого ухудшения разрешения
/// См.: R. B. Blackman, J. W. Tukey, The Measurement of Power Spectra, 1958; а также F. J. Harris, On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform, Proc. IEEE, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class BlackmanHarrisWindow
{
    private const double a0 = 0.35875;
    private const double a1 = 0.48829;
    private const double a2 = 0.14128;
    private const double a3 = 0.01168;

    private const double pi4 = Consts.pi2 * 2;
    private const double pi6 = Consts.pi2 * 3;

    /// <summary>Значение окна в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N) => Value((double)n / N);

    /// <summary>Значение окна в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x) => 
        a0
        - a1 * Math.Cos(Consts.pi2 * x)
        + a2 * Math.Cos(Consts.pi2 * 2 * x)
        - a3 * Math.Cos(Consts.pi2 * 3 * x);

    /// <summary>Значение параметризованного окна Blackman–Harris в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <param name="a0">Коэффициент a0</param>
    /// <param name="a1">Коэффициент a1</param>
    /// <param name="a2">Коэффициент a2</param>
    /// <param name="a3">Коэффициент a3</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N, double a0, double a1, double a2, double a3) => Value((double)n / N, a0, a1, a2, a3);

    /// <summary>Значение параметризованного окна Blackman–Harris в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <param name="a0">Коэффициент a0</param>
    /// <param name="a1">Коэффициент a1</param>
    /// <param name="a2">Коэффициент a2</param>
    /// <param name="a3">Коэффициент a3</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x, double a0, double a1, double a2, double a3) => 
        a0
        - a1 * Math.Cos(Consts.pi2 * x)
        + a2 * Math.Cos(pi4 * x)
        - a3 * Math.Cos(pi6 * x);
}
