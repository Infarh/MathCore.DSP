namespace MathCore.DSP.WindowFunctions;

/// <summary>Синус‑окно (Sine)</summary>
/// <remarks>
/// Простейшее окно на основе синуса, обладающее узким главным лепестком, но сравнительно высокими боковыми лепестками
/// Применяется в задачах, где важна узкая главная полоса и простота реализации; иногда как построительный блок для других окон
/// См.: F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class SineWindow
{
    /// <summary>Значение синус‑окна в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N) => Value((double)n / N);

    /// <summary>Значение синус‑окна в нормированной точке</summary>
    /// <param name="x">Нормированная позиция в интервале [0;1]</param>
    /// <returns>Значение оконной функции в точке x</returns>
    public static double Value(double x) => Math.Sin(Math.PI * x);
}
