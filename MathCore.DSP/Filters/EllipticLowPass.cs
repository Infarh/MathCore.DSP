using static System.Math;

using static MathCore.Polynom.Array;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters;

public class EllipticLowPass : EllipticFilter
{
    public static int GetOrder(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.01)
    {
        var kEps = Sqrt((1 / (Gp * Gp) - 1) / (1 / (Gs * Gs) - 1));
        var pi_dt = PI * dt;
        var kW = Tan(fs * pi_dt) / Tan(fp * pi_dt);

        var N = (int)Ceiling(T(kEps) * K(kW) / (K(kEps) * T(kW)));
        return N;
    }

    public static (IEnumerable<Complex> Zeros, IEnumerable<Complex> Poles) GetPoles(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.01)
    {
        var N = GetOrder(fp, fs, Gp, Gs);

        var EpsP = 1 / (Gp * Gp) - 1;
        var EpsS = 1 / (Gs * Gs) - 1;
        var zeros = EnumNormedZeros(N, EpsP, EpsS);
        var poles = EnumNormedPoles(N, EpsP, EpsS);

        var low_pass_zeros = TransformToLowPass(zeros, fp);
        var low_pass_poles = TransformToLowPass(poles, fp);

        var z_zeros = ToZ(low_pass_zeros, dt);
        if (N.IsOdd())
            z_zeros = z_zeros.AppendFirst(-1);
        var z_poles = ToZ(low_pass_poles, dt);

        return (z_zeros, z_poles);
    }

    /// <summary>Инициализация коэффициентов передаточной функции Эллиптического фильтра</summary>
    /// <returns>Кортеж с коэффициентами полинома числителя и знаменателя передаточной функции</returns>
    private static (double[] A, double[] B) Initialize(Specification opt)
    {
        var k_W = 1 / opt.kW;
        var k_eps = 1 / opt.kEps;

        var N = (int)Ceiling(T(k_eps) * K(k_W) / (K(k_eps) * T(k_W))); // Порядок фильтра

        var eps_p = opt.EpsP;
        var (zeros, poles) = EnumNormedZerosPoles(N, eps_p, opt.EpsS);

        var z_poles = ToZArray(poles, opt.dt, opt.Wp);
        var z_zeros_enum = ToZ(zeros, opt.dt, opt.Wp);
        var is_odd = N.IsOdd();
        if (is_odd) z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = z_zeros_enum.ToArray();

        var k_zeros = z_zeros.Multiply(z => 1 - z);
        var k_poles = z_poles.Multiply(z => 1 - z);
        var k0 = is_odd 
            ? 1 
            : 1 / Sqrt(1 + eps_p * eps_p);

        var g_norm = k0 * (k_poles / k_zeros).Abs;

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