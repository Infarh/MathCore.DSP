using System.Data;

namespace MathCore.DSP.Filters;

public class EllipticBandStop : EllipticFilter
{
    private static Specification GetSpecification(
        double dt,
        double fpl,
        double fsl,
        double fsh,
        double fph,
        double Gp = 0.891250938,
        double Gs = 0.031622777)
    {
        var Fpl = ToAnalogFrequency(fpl, dt);
        var Fsl = ToAnalogFrequency(fsl, dt);
        var Fsh = ToAnalogFrequency(fsh, dt);
        var Fph = ToAnalogFrequency(fph, dt);

        var Wpl = Consts.pi2 * Fpl;
        var Wsl = Consts.pi2 * Fsl;
        var Wsh = Consts.pi2 * Fsh;
        var Wph = Consts.pi2 * Fph;

        var Wc = Wsl * Wsh;
        var dW = Wsh - Wsl;

        var Wp = Wc / Wph > Wpl
            ? Wph
            : Wpl;
        var W0 = Math.Abs(dW * Wp / (Wc - Wp.Pow2()));
        const double W1 = 1;
        var Fp = W0 / Consts.pi2;
        const double Fs = 1 / Consts.pi2;

        var fp = ToDigitalFrequency(Fp, dt);
        var fs = ToDigitalFrequency(Fs, dt);

        return new Specification(dt, fp, fs, Gp, Gs);
    }

    private static (double[] A, double[] B) Initialize(double fsl, double fsh, Specification Spec)
    {
        var kW = 1 / Spec.kW;
        var k_eps = 1 / Spec.kEps;

        var N = (int)Math.Ceiling(T(k_eps) * K(kW) / (K(k_eps) * T(kW))); // Порядок фильтра

        var (zeros, poles) = GetNormedZeros(N, Spec.EpsP, Spec.EpsS, Spec.Wp);

        var dt = Spec.dt;
        var Fsl = ToAnalogFrequency(fsl, dt);
        var Fsh = ToAnalogFrequency(fsh, dt);

        var pzf_zeros = TransformToBandStop(zeros, Fsl, Fsh);
        var pzf_poles = TransformToBandStop(poles, Fsl, Fsh);

        if (N.IsOdd())
        {
            var wc_sqrt = Consts.pi2 * (fsl * fsh).Sqrt(); // странно что не Fsl и Fsh
            pzf_zeros = pzf_zeros.AppendLast(
                (0, +wc_sqrt),
                (0, -wc_sqrt));
        }

        var z_zeros = ToZArray(pzf_zeros, dt);
        var z_poles = ToZArray(pzf_poles, dt);

        var G_norm = (N.IsOdd() ? 1 : Spec.Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    public EllipticBandStop(
        double dt,
        double fpl,
        double fsl,
        double fsh,
        double fph,
        double Gp = 0.891250938,
        double Gs = 0.031622777)
        : this(fsl, fsh, GetSpecification(dt, fpl, fsl, fsh, fph, Gp, Gs))
    {

    }

    private EllipticBandStop(double fsl, double fsh, Specification Spec) : this(Initialize(fsl, fsh, Spec), Spec) { }

    private EllipticBandStop((double[] A, double[] B) Polynoms, Specification Spec) : this(Polynoms.B, Polynoms.A, Spec) { }

    private EllipticBandStop(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
}