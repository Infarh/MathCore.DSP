using static System.Math;

using static MathCore.Consts;

using Signal = System.Func<double, double>;
using IntSpectrum = System.Func<int, MathCore.Complex>;
using DoubleSpectrum = System.Func<double, MathCore.Complex>;

namespace MathCore.DSP.Fourier;

public static class DoubleFuncFT
{
    /// <summary>Выполнить преобразование Фурье</summary>
    /// <param name="s">Вещественная функция</param>
    /// <param name="t1">Начало интервала</param>
    /// <param name="t2">Конец интервала</param>
    /// <param name="IsInverse">Обратное преобразование</param>
    /// <param name="dt">Шаг численного расчёта</param>
    /// <returns>Спектр</returns>
    public static DoubleSpectrum GetFourierTransformation(
        this Signal s,
        double t1,
        double t2,
        bool IsInverse = false,
        double dt = 1e-4)
    {
        var delta_t = t2 - t1;
        var N = (int)Abs(delta_t / dt);
        dt = delta_t / N;
        var w = IsInverse ? -pi2 : pi2;

        return f =>
        {
            var pif = IsInverse ? -w * f : w * f;
            var t = t1;
            var val = s(t);
            var arg = pif * t;
            var (p, q) = (val * Cos(arg), val * Sin(arg));

            var re_s = .0;
            var im_s = .0;
            for (var i = 0; i < N; i++)
            {
                val = s(t);
                arg = pif * t;
                re_s += p + (p = val * Cos(arg));
                im_s += q + (q = val * Sin(arg));
                t += dt;
            }
            return new Complex(re_s * .5 * dt, im_s * .5 * dt);
        };
    }

    /// <summary>Выполнить преобразование Фурье</summary>
    /// <param name="s">Комплексная функция</param>
    /// <param name="t1">Начало интервала</param>
    /// <param name="t2">Конец интервала</param>
    /// <param name="IsInverse">Обратное преобразование</param>
    /// <param name="dt">Шаг численного расчёта</param>
    /// <returns>Спектр</returns>
    public static DoubleSpectrum GetFourierTransformation(
        this DoubleSpectrum s,
        double t1,
        double t2,
        bool IsInverse = false,
        double dt = 1e-4)
    {
        var delta_t = t2 - t1;
        var N = (int)Abs(delta_t / dt);
        dt = delta_t / N;
        var w = IsInverse ? -pi2 : pi2;

        return f =>
        {   //todo: проверить алгоритм!
            var pif = w * f;
            var t = t1;
            var z = s(t).Rotate(pif * t);

            var (re, im) = z;
            for (var i = 0; i < N; i++)
            {
                z = s(t).Rotate(pif * t);
                re += z.Re;
                im += z.Im;
                t += dt;
            }
            return .5 * dt * new Complex(re, im);
        };
    }

    /// <summary>Спектр по целочисленным значением частот</summary>
    /// <param name="s">Вещественная функция</param>
    /// <param name="t1">Начало интервала</param>
    /// <param name="t2">Конец интервала</param>
    /// <param name="IsInverse">Обратное преобразование</param>
    /// <param name="dt">Шаг численного расчёта</param>
    /// <returns>Спектр по целочисленным значениям частот</returns>
    public static IntSpectrum GetFourierSpectrum(
        this Signal s,
        double t1,
        double t2,
        bool IsInverse = false,
        double dt = 1e-4)
    {
        var delta_t = t2 - t1;
        var N = (int)Abs(delta_t / dt);
        var ss = new double[N];
        dt = delta_t / N;
        for (var n = 0; n < N; n++)
            ss[n] = s(t1 + n * dt);
        return ss.GetFourierTransformation(IsInverse);
    }

    /// <summary>Спектр по целочисленным значением частот</summary>
    /// <param name="s">Комплексная функция</param>
    /// <param name="t1">Начало интервала</param>
    /// <param name="t2">Конец интервала</param>
    /// <param name="IsInverse">Обратное преобразование</param>
    /// <param name="dt">Шаг численного расчёта</param>
    /// <returns>Спектр по целочисленным значениям частот</returns>
    public static IntSpectrum GetFourierSpectrum(
        this DoubleSpectrum s,
        double t1,
        double t2,
        bool IsInverse = false,
        double dt = 1e-4)
    {
        var delta_t = t2 - t1;
        var N = (int)Abs(delta_t / dt);

        var ss = new Complex[N];
        dt = delta_t / N;
        for (var n = 0; n < N; n++)
            ss[n] = s(t1 + n * dt);
        return ss.GetFourierTransformation(IsInverse);
    }
}