using System;
using System.Linq;
using System.Runtime.CompilerServices;

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
        /// <param name="dt">Период дискретизации</param>
        /// <param name="fp">Частота пропускания</param>
        /// <param name="fs">Частота заграждения</param>
        /// <param name="Gp">Затухание в полосе пропускания</param>
        /// <param name="Gs">Затухание в полосе заграждения</param>
        /// <returns>Кортеж с коэффициентами полинома числителя и знаменателя передаточной функции</returns>
        private static (double[] A, double[] B) Initialize(double dt, double fp, double fs, double Gp, double Gs)
        {
            if (!(fp < fs)) throw new InvalidOperationException("Частота пропускания должна быть меньше частоты подавления");
            if (!(fp < 1 / (2 * dt))) throw new InvalidOperationException();

            var Rp = -Gp.In_dB();
            var Rs = -Gs.In_dB();

            // Рассчитываем частоты цифрового фильтра
            var Fp = ToAnalogFrequency(fp, dt);
            var Fs = ToAnalogFrequency(fs, dt);

            // Круговые частоты
            var Wp = Consts.pi2 * Fp;

            // Допуск на АЧХ в интервале пропускания
            var eps_p = (Math.Pow(10, Rp / 10) - 1).Sqrt();
            // Допуск на АЧХ в интервале подавления
            var eps_s = (Math.Pow(10, Rs / 10) - 1).Sqrt();

            var k_W = fp / fs;
            var k_eps = eps_p / eps_s;

            var K_w = K(k_W);
            var T_w = T(k_W);
            var K_eps = K(k_eps);
            var T_eps = T(k_eps);

            // Оценка снизу порядка фильтра
            var double_N = T_eps * K_w / (K_eps * T_w);

            var N = (int)Math.Ceiling(double_N); // Порядок фильтра

            var L = N / 2;  // Число комплексно сопряжённых полюсов
            var r = N % 2;  // Число (0 или 1) действительных полюсов - (чётность фильтра)

            // Эллиптический модуль
            var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);

            var m = (1 - k_eps * k_eps).Sqrt();
            var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));
            k_W = (1 - kp * kp).Sqrt();
            var v0_complex = sn_inverse((0, 1 / eps_p), k_eps) / N;

            var zeros = new Complex[N - r]; // Массив нулей (на r меньше числа полюсов)
            var poles = new Complex[N];     // Массив полюсов

            // Если фильтр нечётный, то первым полюсом будет действительный полюс
            if (r != 0) poles[0] = Complex.i * sn_uk(v0_complex, k_W);
            for (var i = 0; i < L; i++)
            {
                // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
                var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_W);

                poles[r + 2 * i] = (-p_re, p_im);
                poles[r + 2 * i + 1] = poles[r + 2 * i].ComplexConjugate;

                var p0_im = 1 / (k_W * cd_uk(u[i], k_W));
                zeros[2 * i] = (0, p0_im);
                zeros[2 * i + 1] = zeros[2 * i].ComplexConjugate;
            }

            //var z_zeros = zeros.ToArray(z => ToZ(z * Wp, dt));
            //var z_poles = poles.ToArray(z => ToZ(z * Wp, dt));

            var z_zeros = ToZArray(zeros, dt, Wp);
            var z_poles = ToZArray(poles, dt, Wp);

            if (r > 0)
            {
                Array.Resize(ref z_zeros, z_zeros.Length + 1);
                Array.Copy(z_zeros, 0, z_zeros, 1, z_zeros.Length - 1);
                z_zeros[0] = -1;
            }

            var G_norm = (r > 0 ? 1 : 1 / (1 + eps_p * eps_p).Sqrt())
                / (z_zeros.Aggregate(Complex.Real, (Z, z) => Z * (1 - z))
                    / z_poles.Aggregate(Complex.Real, (Z, z) => Z * (1 - z))).Abs;

            var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
            var A = GetCoefficientsInverted(z_poles).ToRe();

            return (A, B);
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
            : this(Initialize(dt, fp, fs, Gp, Gs))
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