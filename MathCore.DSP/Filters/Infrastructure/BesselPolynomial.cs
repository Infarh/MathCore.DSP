namespace MathCore.DSP.Filters.Infrastructure;

/// <summary>Класс для вычисления полиномов Бесселя</summary>
internal static class BesselPolynomial
{
    /// <summary>Вычисляет значение полинома Бесселя n-й степени в точке x</summary>
    /// <param name="x">Аргумент полинома</param>
    /// <param name="n">Степень полинома (n ≥ 0)</param>
    /// <returns>Значение полинома Бесселя n-й степени в точке x</returns>
    /// <exception cref="ArgumentOutOfRangeException">Если n &lt; 0</exception>
    public static double Th(double x, int n) => n switch
    {
        < 0 => throw new ArgumentOutOfRangeException(nameof(n), n, "n должно быть >= 0"),
        0 => 1,
        1 => x + 1,
        //2 => x * x + 3 * x + 3,
        //3 => x * x * x + 6 * x * x + 15 * x + 15,
        //4 => x * x * x * x + 10 * x * x * x + 45 * x * x + 105 * x + 105,
        //5 => x * x * x * x * x + 15 * x * x * x * x + 105 * x * x * x + 420 * x * x + 945 * x + 945,
        2 => x * (x + 3) + 3,                                           // x^2 +  3x   +   3
        3 => x * (x * (x + 6) + 15) + 15,                               // x^3 +  6x^2 +  15x   + 15
        4 => x * (x * (x * (x + 10) + 45) + 105) + 105,                 // x^4 + 10x^3 +  45x^2 + 105x   + 105
        5 => x * (x * (x * (x * (x + 15) + 105) + 420) + 945) + 945,    // x^5 + 15x^4 + 105x^3 + 420x^2 + 945x + 945
        _ => (2 * n - 1) * Th(x, n - 1) + x * x * Th(x, n - 2)
    };
}
