namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно Дольфа–Чебышева (Chebyshev, Dolph–Chebyshev)</summary>
/// <remarks>
/// Параметрическое окно с равноволновыми (equiripple) боковыми лепестками, позволяющее задавать их уровень через параметр gamma (в дБ)
/// Применяется в задачах, где требуется строгий контроль уровня боковых лепестков при конкретной ширине главного лепестка; подходит для узкополосного спектрального анализа и синтеза FIR‑фильтров
/// Источник: C. L. Dolph, A current distribution for broadside arrays which optimizes the relationship between beam width and side-lobe level, Proc. IRE, 1946; также F. J. Harris, 1978, https://ieeexplore.ieee.org/document/1455100
/// </remarks>
public static class ChebyshevWindow
{
    /// <summary>Значение окна Чебышева в дискретной точке</summary>
    /// <param name="n">Номер отсчёта (0 ≤ n &lt; N)</param>
    /// <param name="N">Размер окна</param>
    /// <param name="gamma">Уровень боковых лепестков в дБ</param>
    /// <returns>Значение оконной функции в точке n</returns>
    public static double Value(int n, int N, int gamma)
    {
        var q = 10.Pow(gamma / 20);
        
        var (m, d) = N % 2 == 0
            ? (N / 2 - 1, 0.5)  // N - чётное
            : ((N - 1) / 2, 0); // N - нечётное

#if NET8_0_OR_GREATER
        var b = Math.Cosh(Math.Acosh(q) / (N - 1));
#else
        var b = Math.Cosh(MathEx.Hyperbolic.Acosh(q) / (N - 1));
#endif

        var pin = Consts.pi / N;
        var p2in = Consts.pi2 / N;

        var sum = 0d;
        for(var i = 1; i <= m; i++)
        {
            var x1 = b * Math.Cos(pin * i);
            var x2 = Math.Cos(p2in * (n - m - d) * i);
            var c = C(N - 1, x1 * x2);
            sum += c;
        }

        var result = q + 2 * sum;

        return result;

        static double C(int N1, double x) => Math.Cos(N1 * Math.Acos(x));
    }
}
