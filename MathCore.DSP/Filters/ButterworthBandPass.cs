using static System.Math;

using static MathCore.Polynom.Array;

// ReSharper disable InconsistentNaming
// ReSharper disable MemberCanBePrivate.Global
// ReSharper disable UnusedAutoPropertyAccessor.Global

namespace MathCore.DSP.Filters;

public class ButterworthBandPass : ButterworthFilter
{
    private static (double[] A, double[] B) GetPolynoms(Specification Spec, double fmin, double fmax)
    {
        var fd = Spec.fd;
        var Fmin = ToAnalogFrequency(fmin, Spec.dt);
        var Fmax = ToAnalogFrequency(fmax, Spec.dt);

        var N = (int)Ceiling(Log(Spec.kEps) / Log(Spec.kW));
        var poles = GetNormPoles(N, Spec.EpsP);

        //var ppf_zeros = Enumerable.Repeat(new Complex(), N);
        var ppf_poles = TransformToBandPassPoles(poles, Fmin, Fmax).ToArray();

        var ppf_zeros_z = Enumerable.Repeat(Complex.ReValue(1), N).AppendLast(Enumerable.Repeat(Complex.ReValue(-1), N)).ToArray();
        var ppf_poles_z = ToZArray(ppf_poles, Spec.dt);

        var kz0 = Complex.Real;
        var dW = Consts.pi2 * (Fmax - Fmin);
        var fd2 = 2 * Spec.fd;
        for (var i = 0; i < ppf_poles.Length; i += 2)
        {
            kz0 *= fd2 * dW;
            kz0 /= (fd2 - ppf_poles[i]);
            kz0 /= (fd2 - ppf_poles[i + 1]);
        }
        kz0 /= Spec.EpsP;

        var B = GetCoefficientsInverted(ppf_zeros_z).ToArray(v => (v * kz0).Re);
        var A = GetCoefficientsInverted(ppf_poles_z).ToRe()!;

        return (A, B);
    }

    public Interval F { get; }

    public ButterworthBandPass(
        double dt,
        double fp,
        double fs,
        double fmin,
        double fmax,
        double Gp = 0.891250938,
        double Gs = 0.031622777)
        : this(GetSpecification(dt, fp, fs, Gp, Gs), fmin, fmax)
    { }

    private ButterworthBandPass(Specification Spec, double fmin, double fmax)
        : this(GetPolynoms(Spec, fmin, fmax), Spec) => F = (fmin, fmax);

    private ButterworthBandPass((double[] A, double[] B) Polynoms, Specification Spec) : this(Polynoms.B, Polynoms.A, Spec) { }

    private ButterworthBandPass(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
}