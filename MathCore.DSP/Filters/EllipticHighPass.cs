using System.Net.NetworkInformation;

using static System.Array;
using static System.Math;

namespace MathCore.DSP.Filters;

public class EllipticHighPass : EllipticFilter
{
    public static (double[] A, double[] B) GetPolynoms(Specification opt)
    {
        var k_w = opt.kw;
        var k_eps = 1 / opt.kEps;

        var Kw = K(k_w);
        var Tw = T(k_w);
        var K_eps = K(k_eps);
        var T_eps = T(k_eps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N); // Порядок фильтра

        var (zeros, poles) = GetNormedZeros(N, opt.EpsP, opt.EpsS);

        // Масштабируем полюса на требуемую частоту пропускания
        var Fp = opt.Fp;
        var translated_zeros = TransformToHighPass(zeros, Fp);
        var translated_poles = TransformToHighPass(poles, Fp);

        var is_odd = N.IsOdd();
        if (is_odd)
            translated_zeros = translated_zeros.Prepend(0);

        var zz = translated_zeros.ToArray();
        var pp = translated_poles.ToArray();

        // Переходим из p-плоскости в z-плоскость
        var dt = opt.dt;
        var z_zeros = ToZArray(zz, dt);
        var z_poles = ToZArray(pp, dt);
        //var z_poles = translated_poles.ToArray(p => ToZ(p, dt));
        // Вычисляем нормирующий множитель
        //var kz = GetNormalizeCoefficient(translated_poles, dt);

        var Gz0 = (z_zeros.Multiply(z => 1 + z) / z_poles.Multiply(z => 1 + z)).Abs;
        var G_norm = (is_odd ? 1 : opt.Gp) / Gz0;

        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToArray(v => v.Re * G_norm);
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