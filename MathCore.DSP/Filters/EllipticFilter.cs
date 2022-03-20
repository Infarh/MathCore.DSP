using System.Runtime.CompilerServices;
using System.Runtime.Serialization;

using static System.Linq.Enumerable;

using static MathCore.SpecialFunctions.EllipticJacobi;

#nullable enable
namespace MathCore.DSP.Filters;

/// <summary>Эллиптический фильтр</summary>
[KnownType(typeof(EllipticLowPass))]
//[KnownType(typeof(EllipticBandPass))]
public abstract class EllipticFilter : AnalogBasedFilter
{
    protected static (Complex[] Zeros, Complex[] Poles) GetNormedZeros(int N, double EpsP, double EpsS)
    {
        var k_eps = EpsP / EpsS;

        var (L, r) = N.GetDivMod(2);

        // Эллиптический модуль
            var u = Range(1, L).ToArray(i => (2 * i - 1d) / N);

            var m = (1 - k_eps.Pow2()).Sqrt();
            var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Pow2().Pow2());
            var k_w = (1 - kp.Pow2()).Sqrt();
            var v0_complex = sn_inverse((0, 1 / EpsP), k_eps) / N;

            // нулей всегда чётное число (всегда парные)
            var zeros = new Complex[N - r]; // Массив нулей (на r меньше числа полюсов)
            var poles = new Complex[N];     // Массив полюсов

            // Если фильтр нечётный, то первым полюсом будет действительный полюс
            if (r != 0) poles[0] = Complex.i * sn_uk(v0_complex, k_w);
            for (var i = 0; i < L; i++)
            {
                // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
                var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_w);

                poles[r + 2 * i] = (-p_re, p_im);
                poles[r + 2 * i + 1] = poles[r + 2 * i].ComplexConjugate;

                var p0_im = 1 / (k_w * cd_uk(u[i], k_w));
                zeros[2 * i] = (0, p0_im);
                zeros[2 * i + 1] = zeros[2 * i].ComplexConjugate;
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