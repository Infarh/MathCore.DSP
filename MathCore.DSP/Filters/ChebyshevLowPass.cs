using System;
using System.ComponentModel;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

using static MathCore.Polynom.Array;

namespace MathCore.DSP.Filters
{
    /// <summary>Фильтр Чебышева нижних частот</summary>
    public class ChebyshevLowPass : ChebyshevFilter
    {
        private static (double[] A, double[] B) InitializeI(Specification opt)
        {
            var N = (int)Math.Ceiling(arch(opt.kEps) / arch(opt.kW)); // Порядок фильтра

            var poles = GetNormedPolesI(N, opt.EpsP);
            var z_poles = ToZArray(poles, opt.dt, opt.Wp);

            var A = GetCoefficientsInverted(z_poles).ToRe();

            var g_norm = (N.IsOdd() ? 1 : opt.Gp)
                / (2.Power(N) / z_poles.Multiply(z => 1 - z).Re);

            var B = Enumerable
               .Range(0, N + 1)
               .ToArray(i => g_norm * SpecialFunctions.BinomialCoefficient(N, i));

            return (A, B);
        }

        private static (double[] A, double[] B) InitializeII(Specification opt)
        {
            var N = (int)Math.Ceiling(arch(opt.kEps) / arch(opt.kW)); // Порядок фильтра
            var (zeros, poles) = GetNormedPolesII(N, opt.EpsS);

            var z_zeros = N.IsEven()
               ? ToZArray(zeros, opt.dt, opt.Wp)
               : ToZ(zeros, opt.dt, opt.Wp).Prepend(-1).ToArray();
            var z_poles = ToZArray(poles, opt.dt, opt.Wp);

            var B = GetCoefficientsInverted(z_zeros).ToRe();
            var A = GetCoefficientsInverted(z_poles).ToRe();

            var g_norm = 1 / (z_zeros.Multiply(z => 1 - z).Re / z_poles.Multiply(z => 1 - z).Re);

            for (var i = 0; i < B!.Length; i++)
                B[i] *= g_norm;

            return (A, B);
        }

        private static (double[] A, double[] B) InitializeIICorrected(Specification opt)
        {
            var N = (int)Math.Ceiling(arch(opt.kEps) / arch(opt.kW)); // Порядок фильтра
            var (zeros, poles) = GetNormedPolesII(N, opt.EpsS);

            var kw = opt.kw;
            var z_zeros = N.IsEven()
                ? ToZArray(zeros, opt.dt, opt.Wp * kw)
                : ToZ(zeros, opt.dt, opt.Wp * kw).Prepend(-1).ToArray();
            var z_poles = ToZArray(poles, opt.dt, opt.Wp * kw);

            var B = GetCoefficientsInverted(z_zeros).ToRe();
            var A = GetCoefficientsInverted(z_poles).ToRe();

            var g_norm = 1 / (z_zeros.Multiply(z => 1 - z).Re / z_poles.Multiply(z => 1 - z).Re);

            for (var i = 0; i < B!.Length; i++)
                B[i] *= g_norm;

            return (A, B);
        }

        public ChebyshevType FilterType { get; }

        /// <summary>Инициализация нового фильтра Чебышева нижних частот</summary>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="fp">Частота пропускания</param>
        /// <param name="fs">Частота заграждения</param>
        /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
        /// <param name="Gs">Затухание в полосе заграждения (0.005623413 = -45 дБ)</param>
        /// <param name="Type">Тип (род) фильтра чебышева</param>
        public ChebyshevLowPass(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.005623413, ChebyshevType Type = ChebyshevType.I)
            : this(GetSpecification(dt, fp, fs, Gp, Gs), Type) => FilterType = Type;

        private ChebyshevLowPass(Specification opt, ChebyshevType Type)
            : this(Type switch
            {
                ChebyshevType.I => InitializeI(opt),
                ChebyshevType.II => InitializeII(opt),
                ChebyshevType.IICorrected => InitializeIICorrected(opt),
                _ => throw new InvalidEnumArgumentException(nameof(Type), (int)Type, typeof(ChebyshevType))
            }) { }

        private ChebyshevLowPass((double[] A, double[] B) config) : base(config.B, config.A) { }
    }
}
