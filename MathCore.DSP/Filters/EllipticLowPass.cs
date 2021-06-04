using System;
using System.Linq;
using System.Runtime.CompilerServices;

using static System.Array;
using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions.EllipticJacobi;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters
{
    public class EllipticLowPass : EllipticFilter
    {
        /// <summary>Полный эллиптический интеграл</summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double K(double k) => FullEllipticIntegral(k);

        /// <summary>Полный комплиментарный эллиптический интеграл</summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double T(double k) => FullEllipticIntegralComplimentary(k);

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

            var G_norm = (is_odd ? 1 : 1 / (1 + opt.EpsP.Pow2()).Sqrt())
                / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

            var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
            var A = GetCoefficientsInverted(z_poles).ToRe();

            return (A!, B!);
        }

        public double fp { get; }
        public double fs { get; }
        public double dt { get; }
        public double Gp { get; }
        public double Gs { get; }

        /// <summary>Инициализация нового Эллиптического фильтра нижних частот</summary>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="fp">Частота пропускания</param>
        /// <param name="fs">Частота заграждения</param>
        /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
        /// <param name="Gs">Затухание в полосе заграждения (0.005623413 = -45 дБ)</param>
        public EllipticLowPass(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.005623413)
            : this(Initialize(GetSpecification(dt, fp, fs, Gp, Gs)))
        {
            this.fp = fp;
            this.fs = fs;
            this.dt = dt;
            this.Gp = Gp;
            this.Gs = Gs;
        }

        /// <summary>Инициализация нового Эллиптического фильтра</summary>
        /// <param name="config">Кортеж, содержащий массив коэффициентов полинома числителя и знаменателя</param>
        private EllipticLowPass((double[] A, double[] B) config) : base(config.B, config.A) { }
    }
}