using System;
using System.Linq;

namespace MathCore.DSP.Filters
{
    /// <summary>Низкочастотный фильтр Баттерворта</summary>
    public class ButterworthLowPass : ButterworthFilter
    {
        /// <summary>Конфигурация низкочастотного фильтра Баттерворта, содержащая набор коэффициентов прямой и обратной связи</summary>
        private struct ButterworthLowPassConfiguration
        {
            /// <summary>Коэффициенты прямой связи фильтра</summary>
            public readonly double[] A;

            /// <summary>Коэффициенты обратной связи фильтра</summary>
            public readonly double[] B;

            /// <summary>Инициализация новой конфигурации низкочастотного фильтра Баттерворта</summary>
            /// <param name="fp">Граничная частота пропускания</param>
            /// <param name="fs">Граничная частота подавления</param>
            /// <param name="dt">Период дискретизации</param>
            /// <param name="Gp">Допустимое затухание в полосе пропускания</param>
            /// <param name="Gs">Допустимое затухание в полосе подавления</param>
            public ButterworthLowPassConfiguration(double fp, double fs, double dt, double Gp, double Gs)
            {
                if(!(fp < fs)) throw new InvalidOperationException("Частота пропускания должна быть меньше частоты подавления");
                if(!(fp < 1/ (2 * dt))) throw new InvalidOperationException();

                var Rp = -Gp.In_dB();
                var Rs = -Gs.In_dB();

                var eps_p = Math.Sqrt(Math.Pow(10, Rp / 10) - 1);
                var eps_s = Math.Sqrt(Math.Pow(10, Rs / 10) - 1);


                var Fp = ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
                var Fs = ToAnalogFrequency(fs, dt);  // Частота подавления аналогового пропита

                var Wp = Consts.pi2 * Fp;

                var k_eps = eps_s / eps_p;
                var k_W = Fs / Fp;

                var N = (int)Math.Ceiling(Math.Log(k_eps) / Math.Log(k_W));
                var r = N % 2;

                var alpha = Math.Pow(eps_p, -1d / N);

                var th0 = Math.PI / N;

                var poles = new Complex[N];
                if (r != 0) poles[0] = -alpha;
                for (var i = r; i < poles.Length; i += 2)
                {
                    var w = th0 * (i + 1 - r - 0.5);
                    var sin = -alpha * Math.Sin(w);
                    var cos = alpha * Math.Cos(w);
                    poles[i] = new Complex(sin, cos);
                    poles[i + 1] = new Complex(sin, -cos);
                }

                var translated_poles = poles.ToArray(p => p * Wp);
                var z_poles = translated_poles.ToArray(p => ToZ(p, dt));
                var kz = GetNomalizeCoefficient(translated_poles, dt);
                var WpN = Math.Pow(Wp, N);
                var k = WpN * kz / eps_p;
                B = new double[N + 1];
                for (var i = 0; i < B.Length; i++)
                    B[i] = k * SpecialFunctions.BinomialCoefficient(N, i);
                A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();
            }
        }

        /// <summary>Инициализация нового фильтра Баттерворта нижних частот</summary>
        /// <param name="fp">Частота пропускания</param>
        /// <param name="fs">Частота заграждения</param>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="Gp">Затухание в полосе пропускания</param>
        /// <param name="Gs">Затухание в полосе заграждения</param>
        public ButterworthLowPass(double fp, double fs, double dt, double Gp = 0.891250938, double Gs = 0.031622777)
            : this(new ButterworthLowPassConfiguration(fp, fs, dt, Gp, Gs)) { }

        private ButterworthLowPass(ButterworthLowPassConfiguration config) : base(config.B, config.A) { }
    }
}
