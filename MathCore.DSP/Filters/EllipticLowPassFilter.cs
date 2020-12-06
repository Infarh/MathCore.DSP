using System;
using System.Linq;

namespace MathCore.DSP.Filters
{
    public class EllipticLowPassFilter : EllipticFilter
    {
        /// <summary>Полный эллиптический интеграл</summary>
        private static double K(double k) => SpecialFunctions.EllipticJacobi.FullEllipticIntegral(k);

        /// <summary>Полный комплиментарный эллиптический интеграл</summary>
        private static double T(double k) => SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(k);

        /// <summary>Инициализация коэффициентов передаточной функции Эллиптического фильтра</summary>
        /// <param name="fp">Частота пропускания</param>
        /// <param name="fs">Частота заграждения</param>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="Gp">Затухание в полосе пропускания</param>
        /// <param name="Gs">Затухание в полосе заграждения</param>
        /// <returns>Кортеж с коэффициентами полинома числителя и знаменателя передаточной функции</returns>
        private static (double[] A, double[] B) Initialize(double fp, double fs, double dt, double Gp, double Gs)
        {
            if (!(fp < fs)) throw new InvalidOperationException("Частота пропускания должна быть меньше частоты подавления");
            if (!(fp < 1 / (2 * dt))) throw new InvalidOperationException();

            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);

            const double Rp = 1;
            const double Rs = 45;

            var eps_p = Math.Sqrt(Math.Pow(10, Rp / 10) - 1);
            var eps_s = Math.Sqrt(Math.Pow(10, Rs / 10) - 1);

            //var k_eps = eps_s / eps_p;
            //var k_W = Fs / Fp;
            //Assert.That.Value(k_eps).IsEqual(349.46669702542425);
            //Assert.That.Value(k_W).IsEqual(1.705275881518411, 2.23e-16);

            var k_W = fp / fs;
            var k_eps = eps_p / eps_s;

            var K_w = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(k_W);
            var T_w = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(k_W);
            var K_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(k_eps);
            var T_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(k_eps);

            // Оценка снизу порядка фильтра
            var double_N = T_eps * K_w / K_eps / T_w;

            var N = (int)Math.Ceiling(double_N); // Порядок фильтра

            var L = N / 2;
            var r = N % 2;

            // Эллиптический модуль
            double U(int i) => (2 * i - 1d) / N;
            var u = new double[L];
            for (var i = 0; i < L; i++)
                u[i] = U(i + 1);

            var m = (1 - k_eps * k_eps).Sqrt();

            var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * SpecialFunctions.EllipticJacobi.sn_uk(ui, m).Power(4));

            k_W = (1 - kp * kp).Sqrt();

            var im_pz = new double[L];
            for (var i = 0; i < L; i++)
                im_pz[i] = 1 / (k_W * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k_W));

            var v0_complex = SpecialFunctions.EllipticJacobi.sn_inverse(new Complex(0, 1 / eps_p), k_eps) / N;

            var Pp = new Complex[N];
            var P0 = new Complex[N - r];

            if (r != 0) Pp[0] = Complex.i * SpecialFunctions.EllipticJacobi.sn_uk(v0_complex, k_W);
            for (var i = 0; i < L; i++)
            {
                var (p_im, p_re) = SpecialFunctions.EllipticJacobi.cd_uk(u[i] - v0_complex, k_W);

                Pp[r + 2 * i] = new Complex(-p_re, p_im);
                Pp[r + 2 * i + 1] = new Complex(-p_re, -p_im);

                var p0_im = 1 / (k_W * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k_W));
                P0[2 * i] = new Complex(0, p0_im);
                P0[2 * i + 1] = new Complex(0, -p0_im);
            }


            var B = Polynom.Array.GetCoefficientsInverted(P0).ToRe();
            var A = Polynom.Array.GetCoefficientsInverted(Pp).ToRe();

            var norm_k = B![B.Length - 1] / A![A.Length - 1];

            for (int i = 0, count = B.Length; i < count; i++)
                B[i] /= norm_k;

            return (A, B);
        }

        public double fp { get; }
        public double fs { get; }
        public double dt { get; }
        public double Gp { get; }
        public double Gs { get; }

        /// <summary>Инициализация нового Эллиптического фильтра нижних частот</summary>
        /// <param name="fp">Частота пропускания</param>
        /// <param name="fs">Частота заграждения</param>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
        /// <param name="Gs">Затухание в полосе заграждения (0.005623413 = -45 дБ)</param>
        public EllipticLowPassFilter(double fp, double fs, double dt, double Gp = 0.891250938, double Gs = 0.005623413)
            : this(Initialize(fp, fs, dt, Gp, Gs))
        {
            this.fp = fp;
            this.fs = fs;
            this.dt = dt;
            this.Gp = Gp;
            this.Gs = Gs;
        }

        /// <summary>Инициализация нового Эллиптического фильтра</summary>
        /// <param name="config">Кортеж, содержащий массив коэффициентов полинома числителя и знаменателя</param>
        private EllipticLowPassFilter((double[] A, double[] B) config) : base(config.B, config.A) { }
    }
}