using static System.Array;
using static System.Math;

namespace MathCore.DSP.Filters;

public class EllipticHighPass : EllipticFilter
{
    public static (double[] A, double[] B) GetPolynoms(Specification opt)
    {
        var k_w = 1 / opt.kw;
        var k_eps = 1 / opt.kEps;

        var N = (int)Ceiling(T(k_eps) * K(k_w) / (K(k_eps) * T(k_w))); // Порядок фильтра

        var (zeros, poles) = GetNormedZeros(N, opt.EpsP, opt.EpsS);

        // Масштабируем полюса на требуемую частоту пропускания
        var Wp = opt.Wp;
        var translated_zeros = TransformToHighPass(zeros, Wp);
        var translated_poles = TransformToHighPass(poles, Wp);

        if (N.IsOdd())
            translated_zeros = translated_zeros.Prepend(1);

        // Переходим из p-плоскости в z-плоскость
        var dt = opt.dt;
        var z_zeros = ToZArray(translated_zeros, dt);
        var z_poles = ToZArray(translated_poles, dt);
        //var z_poles = translated_poles.ToArray(p => ToZ(p, dt));
        // Вычисляем нормирующий множитель
        //var kz = GetNormalizeCoefficient(translated_poles, dt);
        var kz = poles.Multiply();

        var WpN = Wp.Pow(N);
        var k = WpN * kz / opt.EpsP;

        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToArray(v => (v * k).Re);
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    public EllipticHighPass(
        double dt,
        double fs,
        double fp,
        double Gp = 0.891250938,
        double Gs = 0.031622777) 
        : this(GetSpecification(dt, fp, fs, Gp, Gs)) { }

    public EllipticHighPass(Specification Spec) : this(GetPolynoms(Spec), Spec) { }

    public EllipticHighPass((double[] A, double[] B) config, Specification Spec) : base(config.B, config.A, Spec) { }
}