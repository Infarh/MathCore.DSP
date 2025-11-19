namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно Ланцоша (Lanczos)</summary>
/// <remarks>
/// Окно на основе функции sinc, обеспечивающее хорошее подавление боковых лепестков и узкий главный лепесток, часто применяется при интерполяции и ресэмплинге
/// Рекомендуется для задач реконструкции сигналов и изображений, а также при спектральном анализе с акцентом на разрешение
/// См.: C. Lanczos, Evaluation of Noisy Data, J. Soc. Indust. Appl. Math., 1964; дополнительно F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class LanczosWindow
{
    /// <summary>Значение окна Ланцоша в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N) => Value((double)n / N);

    /// <summary>Значение окна Ланцоша в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x) => MathEx.Sinc(Consts.pi2 * x - Consts.pi);

    //public static double Value(int n, int N) => MathEx.Sinc(Consts.pi2 * n / N - Consts.pi);
}
