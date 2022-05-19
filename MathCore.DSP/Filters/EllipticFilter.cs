using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.Serialization;

using static System.Linq.Enumerable;

using static MathCore.SpecialFunctions.EllipticJacobi;
// ReSharper disable InconsistentNaming

#nullable enable
namespace MathCore.DSP.Filters;

/// <summary>Эллиптический фильтр</summary>
[KnownType(typeof(EllipticLowPass))]
[KnownType(typeof(EllipticHighPass))]
[KnownType(typeof(EllipticBandPass))]
[KnownType(typeof(EllipticBandStop))]
public abstract class EllipticFilter : AnalogBasedFilter
{
    protected static IEnumerable<Complex> EnumNormedZeros(int N, double EpsP, double EpsS, double W0 = 1)
    {
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var kEps = EpsP / EpsS;

        var L = N / 2;

        // Эллиптический модуль
        var u = Range(1, L).ToArray(i => (2 * i - 1d) / N);

        var m = (1 - kEps.Pow2()).Sqrt();
        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Pow2().Pow2());
        var k_w = (1 - kp.Pow2()).Sqrt();

        for (var i = 0; i < L; i++)
        {
            var im_value = W0 / (k_w * cd_uk(u[i], k_w));
            yield return new(0, +im_value);
            yield return new(0, -im_value);
        }
    }

    protected static IEnumerable<Complex> EnumNormedPoles(int N, double EpsP, double EpsS, double W0 = 1)
    {
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var kEps = EpsP / EpsS;

        var (L, r) = N.GetDivMod(2);

        var u = Range(1, L).ToArray(i => (2 * i - 1d) / N);

        var m = (1 - kEps.Pow2()).Sqrt();
        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Pow2().Pow2());
        var k_w = (1 - kp.Pow2()).Sqrt();
        var v0_complex = sn_inverse((0, 1 / EpsP), kEps) / N;

        if (r != 0)
            yield return Complex.ImValue(W0) * sn_uk(v0_complex, k_w);
        for (var i = 0; i < L; i++)
        {
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_w) * W0;

            yield return new(-p_re, +p_im);
            yield return new(-p_re, -p_im);
        }
    }

    protected static (IEnumerable<Complex> Zeros, IEnumerable<Complex> Poles) EnumNormedZerosPoles(int N, double EpsP, double EpsS, double W0 = 1)
    {
        var zeros = EnumNormedZeros(N, EpsP, EpsS, W0);
        var poles = EnumNormedPoles(N, EpsP, EpsS, W0);
        return (zeros, poles);
    }

    protected static (Complex[] Zeros, Complex[] Poles) GetNormedZerosPoles(int N, double EpsP, double EpsS, double W0 = 1)
    {
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var kEps = EpsP / EpsS;

        var (L, r) = N.GetDivMod(2);

        // Эллиптический модуль
        var u = Range(1, L).ToArray(i => (2 * i - 1d) / N);

        var m = (1 - kEps.Pow2()).Sqrt();
        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Pow2().Pow2());
        var k_w = (1 - kp.Pow2()).Sqrt();
        var v0_complex = sn_inverse((0, 1 / EpsP), kEps) / N;

        // нулей всегда чётное число (всегда парные)
        var zeros = new Complex[N - r]; // Массив нулей (на r меньше числа полюсов)
        var poles = new Complex[N];     // Массив полюсов

        // Если фильтр нечётный, то первым полюсом будет действительный полюс
        if (r != 0) poles[0] = Complex.ImValue(W0) * sn_uk(v0_complex, k_w);
        for (var i = 0; i < L; i++)
        {
            // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_w) * W0;

            (poles[2 * i + r], poles[2 * i + r + 1]) = Complex.Conjugate(-p_re, p_im);
            (zeros[2 * i], zeros[2 * i + 1]) = Complex.Conjugate(0, W0 / (k_w * cd_uk(u[i], k_w)));
        }

        return (zeros, poles);
    }

    /// <summary>Инициализация нового эллиптического фильтра</summary>
    /// <param name="B">Коэффициенты полинома числителя</param>
    /// <param name="A">Коэффициенты полинома знаменателя</param>
    /// <param name="Spec">Спецификация фильтра</param>
    protected EllipticFilter(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }

    /// <summary>Полный эллиптический интеграл</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    protected static double K(double k) => FullEllipticIntegral(k);
     
    /// <summary>Полный комплиментарный эллиптический интеграл</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    protected static double T(double k) => FullEllipticIntegralComplimentary(k);
}