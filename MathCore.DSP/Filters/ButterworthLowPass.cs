﻿using System;
using System.Linq;

namespace MathCore.DSP.Filters
{
    public abstract class AnalogBasedFilter : IIR
    {
        protected AnalogBasedFilter(double[] B, double[] A) : base(B, A) { }
    }

    public abstract class ButterworthFilter : AnalogBasedFilter
    {
        protected ButterworthFilter(double[] B, double[] A) : base(B, A) { }
    }

    public class ButterworthLowPass : ButterworthFilter
    {
        private struct ButterworthLowPassConfiguration
        {
            public readonly double[] A;
            public readonly double[] B;

            public ButterworthLowPassConfiguration(double fp, double fs, double dt, double Gp, double Gs)
            {
                if(!(fp < fs)) throw new InvalidOperationException("Частота пропускания должна быть меньше частоты подавления");
                if(!(fp < 1/ (2 * dt))) throw new InvalidOperationException();

                //double wp = Consts.pi2 * fp * dt; // 0.628318 рад/с
                //double ws = Consts.pi2 * fs * dt; // 1.884955 рад/с

                var Rp = -Gp.In_dB();
                var Rs = -Gs.In_dB();

                var eps_p = Math.Sqrt(Math.Pow(10, Rp / 10) - 1);
                var eps_s = Math.Sqrt(Math.Pow(10, Rs / 10) - 1);


                var Fp = ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
                var Fs = ToAnalogFrequency(fs, dt);  // Частота подавления аналогового протипа

                var Wp = Consts.pi2 * Fp;
                //var Ws = Consts.pi2 * Fs;

                var k_eps = eps_s / eps_p;
                var k_W = Fs / Fp;

                var N = (int)Math.Ceiling(Math.Log(k_eps) / Math.Log(k_W));
                //var L = N / 2;
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

                var translated_poles = poles.Select(p => p * Wp).ToArray();
                var z_poles = translated_poles.Select(p => ToZ(p, dt)).ToArray();
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
