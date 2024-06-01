using static MathCore.Complex;
using static MathCore.Consts;

namespace MathCore.DSP.Fourier;

/// <summary>Класс преобразования Фурье</summary>
public static class FT
{
    private static (double Sin, double Cos)[] GetCoefficients(int N, bool IsInverse = false)
    {
        var w   = new (double Sin, double Cos)[N];
        var w0  = IsInverse ? -Consts.pi2 / N : Consts.pi2 / N;
        var arg = 0.0;
        for (var i = 0; i < N; i++)
        {
            w[i] =  (Math.Sin(arg), Math.Cos(arg));
            arg  += w0;
        }
        return w;
    }

    /* -------------------------------------------------------------------------------------------- */

    /// <summary>Коэффициент прямого преобразования Фурье</summary>
    /// <param name="N">Число спектральных составляющих</param>
    /// <returns>Комплексное значение коэффициента Фурье</returns>
    public static Complex W(int N) => Exp(1d / N, -pi2 / N);

    /// <summary>Коэффициент обратного преобразования Фурье</summary>
    /// <param name="N">Число спектральных составляющих</param>
    /// <returns>Комплексное значение обратного коэффициента Фурье</returns>
    public static Complex Winv(int N) => Exp(pi2 / N);

    /// <summary>Коэффициент прямого преобразования Фурье для спектральной составляющей</summary>
    /// <param name="k">Номер спектральной составляющей</param>
    /// <param name="N">Число спектральных составляющих</param>
    /// <returns>Комплексное значение коэффициента Фурье</returns>
    public static Complex W(int k, int N) => Exp(1d / N, -pi2 * k / N);

    /// <summary>Коэффициент обратного преобразования Фурье для спектральной составляющей</summary>
    /// <param name="k">Номер спектральной составляющей</param>
    /// <param name="N">Число спектральных составляющих</param>
    /// <returns>Комплексное значение обратного коэффициента Фурье</returns>
    public static Complex Winv(int k, int N) => Exp(1d / N, -pi2 * k / N);

    /* -------------------------------------------------------------------------------------------- */

    public static Complex[] GetFourierTransform(this IEnumerable<double> Values, bool IsInverse = false)
        => Values.ToArray().GetFourierTransform(IsInverse);

    /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
    /// <param name="Values">Массив отсчётов функции</param>
    /// <param name="IsInverse">Обратное преобразование</param>
    public static Complex[] GetFourierTransform(this double[] Values, bool IsInverse = false)
    {
        if (Values is null) throw new ArgumentNullException(nameof(Values));

        var N = Values.Length;
        var spectrum = new Complex[N];
        var w = GetCoefficients(N, IsInverse);

        for (var m = 0; m < N; m++)
        {
            var p = 0.0;
            var q = 0.0;
            for (var n = 0; n < N; n++)
            {
                var v          = Values[n];
                var (sin, cos) =  w[m * n % N];
                p              += v * cos;
                q              += v * sin;
            }
            spectrum[m] = new(p / N, q / N);
        }

        return spectrum;
    }

    /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
    /// <param name="Values">Массив отсчётов функции</param>
    /// <param name="IsInverse">Обратное преобразование</param>
    /// <param name="progress">Метод информирования о прогрессе операции</param>
    public static Complex[] GetFourierTransform(this double[] Values, bool IsInverse, Action<double>? progress)
    {
        if (Values is null) throw new ArgumentNullException(nameof(Values));

        var N = Values.Length;
        var spectrum = new Complex[N];
        var w = GetCoefficients(N, IsInverse);

        for (var m = 0; m < N; m++)
        {
            var p = 0.0;
            var q = 0.0;
            for (var n = 0; n < N; n++)
            {
                var v          = Values[n];
                var (sin, cos) =  w[m * n % N];
                p              += v * cos;
                q              += v * sin;
            }
            spectrum[m] = new(p / N, q / N);
            progress?.Invoke((double)m / N);
        }
        progress?.Invoke(1);

        return spectrum;
    }

    /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
    /// <param name="Values">Массив отсчётов функции</param>
    /// <param name="Inverse">Обратное преобразование</param>
    public static Complex[] GetFourierTransform(this Complex[] Values, bool Inverse = false)
    {
        if (Values is null) throw new ArgumentNullException(nameof(Values));

        var spectrum = new Complex[Values.Length];
        var N = spectrum.Length;
        var w = GetCoefficients(N, Inverse);

        for (var m = 0; m < N; m++)
        {
            var P = 0.0;
            var Q = 0.0;
            for (var n = 0; n < N; n++)
            {
                var (re, im)   =  Values[n];
                var (sin, cos) =  w[n * m % N];
                P              += re * cos - im * sin;
                Q              += im * cos + re * sin;
            }
            spectrum[m] = Inverse ? new(P, Q) : new(P / N, Q / N);
        }
        return spectrum;
    }

    /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
    /// <param name="Values">Массив отсчётов функции</param>
    /// <param name="Inverse">Обратное преобразование</param>
    /// <param name="progress">Метод индикации прогресса выполнения</param>
    public static Complex[] GetFourierTransform(this Complex[] Values, bool Inverse, Action<double>? progress)
    {
        if (Values is null) throw new ArgumentNullException(nameof(Values));

        var spectrum = new Complex[Values.Length];
        var N = spectrum.Length;
        var w = GetCoefficients(N, Inverse);

        for (var m = 0; m < N; m++)
        {
            var P = 0.0;
            var Q = 0.0;
            for (var n = 0; n < N; n++)
            {
                // ReSharper disable once UseDeconstruction
                var (re, im)   =  Values[n];
                var (sin, cos) =  w[n * m % N];
                P              += re * cos - im * sin;
                Q              += im * cos + re * sin;
            }

            progress?.Invoke((double)m / N);
            spectrum[m] = Inverse ? new(P, Q) : new(P / N, Q / N);
        }
        progress?.Invoke(1);
        return spectrum;
    }
}