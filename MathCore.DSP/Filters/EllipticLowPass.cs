using static System.Array;
using static System.Math;

using static MathCore.Polynom.Array;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters;

public class EllipticLowPass : EllipticFilter
{
    /// <summary>Инициализация коэффициентов передаточной функции Эллиптического фильтра</summary>
    /// <returns>Кортеж с коэффициентами полинома числителя и знаменателя передаточной функции</returns>
    private static (double[] A, double[] B) Initialize(Specification opt)
    {
        var k_W = 1 / opt.kw;
        var k_eps = 1 / opt.kEps;

        var N = (int)Ceiling(T(k_eps) * K(k_W) / (K(k_eps) * T(k_W))); // Порядок фильтра

        var (zeros, poles) = GetNormedZeros(N, opt.EpsP, opt.EpsS);

        var z_zeros = ToZArray(zeros, opt.dt, opt.Wp);
        var z_poles = ToZArray(poles, opt.dt, opt.Wp);

        var is_odd = N.IsOdd();
        if (is_odd)
        {
            Resize(ref z_zeros, z_zeros.Length + 1);
            Copy(z_zeros, 0, z_zeros, 1, z_zeros.Length - 1);
            z_zeros[0] = -1;
        }

        var g_norm = (is_odd ? 1 : 1 / (1 + opt.EpsP.Pow2()).Sqrt())
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    /// <summary>Инициализация нового Эллиптического фильтра нижних частот</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="fp">Частота пропускания</param>
    /// <param name="fs">Частота заграждения</param>
    /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
    /// <param name="Gs">Затухание в полосе заграждения (0.01 = -40 дБ)</param>
    public EllipticLowPass(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.01)
        : this(GetSpecification(dt, fp, fs, Gp, Gs)) { }

    public EllipticLowPass(Specification Spec) : this(Initialize(Spec), Spec) { }

    /// <summary>Инициализация нового Эллиптического фильтра</summary>
    /// <param name="Polynoms">Кортеж, содержащий массив коэффициентов полинома числителя и знаменателя</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private EllipticLowPass((double[] A, double[] B) Polynoms, Specification Spec) : base(Polynoms.B, Polynoms.A, Spec) { }
}