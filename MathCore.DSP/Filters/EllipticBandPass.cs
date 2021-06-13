using System;
using System.Linq;
using System.Net.Http.Headers;

using MathCore.Values;

using static System.Array;
using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions.EllipticJacobi;

namespace MathCore.DSP.Filters
{
    public class EllipticBandPass : EllipticFilter
    {
        private static Specification GetSpecification(
            double dt,
            double fsl,
            double fpl,
            double fph,
            double fsh,
            double Gp = 0.891250938,
            double Gs = 0.031622777)
        {
            var fs = fpl * fph > fsl * fsh ? fsh : fsl;
            return new(dt, 1 / Consts.pi2, Abs((fpl * fph - fs * fs) / (fph - fpl) / fs) / Consts.pi2, Gp, Gs);
        }

        private static (double[] A, double[] B) Initialize(double fpl, double fph, Specification Spec)
        {
            var k_W = 1 / Spec.kw;
            var k_eps = 1 / Spec.kEps;

            var N = (int)Ceiling(T(k_eps) * K(k_W) / (K(k_eps) * T(k_W))); // Порядок фильтра

            var (zeros, poles) = GetNormedZeros(N, Spec.EpsP, Spec.EpsS);

            var Fpl = ToAnalogFrequency(fpl, Spec.dt);
            var Fph = ToAnalogFrequency(fph, Spec.dt);

            var ppf_zeros = TransformToBandPassPoles(zeros, Fpl, Fph)
               .Concat(Enumerable.Repeat(new Complex(), poles.Length - zeros.Length))
               .ToArray();
            var ppf_poles = TransformToBandPassPoles(poles, Fpl, Fph).ToArray();

            var is_odd = N.IsOdd();
            var z_zeros = is_odd
                ? ToZ(ppf_zeros, Spec.dt, Spec.Wp).AppendFirst(-1).ToArray()
                : ToZArray(ppf_zeros, Spec.dt, Spec.Wp);
            var z_poles = ToZArray(ppf_poles, Spec.dt, Spec.Wp);

            var z0 = Complex.Exp(-Consts.pi2 * (fpl * fph).Sqrt() / Spec.fd);

            var norm = z_zeros.Aggregate(1d, (Z, z) => Z * (1 - z * z0).Abs)
                / z_poles.Aggregate(1d, (Z, z) => Z * (1 - z * z0).Abs);

            var B = GetCoefficientsInverted(z_zeros).ToArray(b => b.Re / norm);
            var A = GetCoefficientsInverted(z_poles).ToRe();

            return (A!, B!);
        }

        public EllipticBandPass(
            double dt,
            double fsl,
            double fpl,
            double fph,
            double fsh,
            double Gp = 0.891250938,
            double Gs = 0.031622777)
            : this(fsl, fpl, fph, fsh, GetSpecification(dt, fsl, fpl, fph, fsh, Gp, Gs))
        {

        }

        private EllipticBandPass(
            double fsl,
            double fpl,
            double fph,
            double fsh,
            Specification Spec)
            : this(Initialize(fpl, fph, Spec), Spec)
        { }

        private EllipticBandPass((double[] A, double[] B) Polynoms, Specification Spec) : this(Polynoms.B, Polynoms.A, Spec) { }

        private EllipticBandPass(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
    }
}