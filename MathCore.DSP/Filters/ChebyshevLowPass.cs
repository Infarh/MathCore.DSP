﻿using System;
using System.ComponentModel;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

using static MathCore.Polynom.Array;

namespace MathCore.DSP.Filters
{
    public class ChebyshevLowPass : ChebyshevFilter
    {
        private static double arsh(double x) => Math.Log(x + Math.Sqrt(x * x + 1));
        private static double arch(double x) => Math.Log(x + Math.Sqrt(x * x - 1));

        [SuppressMessage("ReSharper", "InconsistentNaming")]
        public enum ChebyshevType : byte
        {
            /// <summary>Фильтр Чебышева первого рода - основной фильтр, пропускающий нижнюю полосу частот</summary>
            I,
            /// <summary>Фильтр Чебышева второго рода, подавляющий верхнюю область частот (выше fp)</summary>
            II,
            /// <summary>Фильтр Чебышева второго рода c коррекцией частотного диапазона, подавляющий верхнюю область частот (выше fs)</summary>
            IICorrected
        }

        private static Complex[] GetAnalogPolesI(int N, double EpsP)
        {
            var r = N % 2;                              // Нечётность порядка фильтра
            var dth = Math.PI / N;                      // Угловой шаг между полюсами
            var beta = arsh(1 / EpsP) / N;
            var sh = Math.Sinh(beta);
            var ch = Math.Cosh(beta);
            var poles = new Complex[N];                 // Массив полюсов фильтра
            if (r != 0) poles[0] = -sh;    // Если порядок фильтра нечётный, то первым добавляем центральный полюс
            for (var i = r; i < poles.Length; i += 2)   // Расчёт полюсов
            {
                var th = dth * (i + 1 - r - 0.5);

                var sin = Math.Sin(th);
                var cos = Math.Cos(th);
                poles[i] = new Complex(-sh * sin, ch * cos);
                poles[i + 1] = poles[i].ComplexConjugate;
            }
            
            return poles;
        }

        private static (Complex[] Zeros, Complex[] Poles) GetAnalogPolesII(int N, double EpsP)
        {
            var r = N % 2;                              // Нечётность порядка фильтра
            var L = (N - r) / 2;                        // Число пар нулей
            var dth = Math.PI / N;                      // Угловой шаг между полюсами
            var beta = arsh(EpsP) / N;
            var shb = Math.Sinh(beta);
            var chb = Math.Cosh(beta);
           
            var poles = new Complex[N];                 // Массив полюсов фильтра
            if (r != 0) poles[0] = -1 / shb;            // Если порядок фильтра нечётный, то первым добавляем центральный полюс
            for (var i = r; i < poles.Length; i += 2)   // Расчёт полюсов
            {
                var th = dth * (i + 1 - r - 0.5);
                var sin = Math.Sin(th);
                var cos = Math.Cos(th);
                var norm = 1 / (sin * sin * shb * shb + cos * cos * chb * chb);
                poles[i] = new Complex(-shb * sin * norm, chb * cos * norm);
                poles[i + 1] = poles[i].ComplexConjugate;
            }

            var zeros = new Complex[L * 2];
            for (var n = 1; n <= L; n++)
            {
                var th = dth * (n - 0.5);
                zeros[2 * n - 2] = new Complex(0, 1 / Math.Cos(th));
                zeros[2 * n - 1] = zeros[2 * n - 2].ComplexConjugate;
            }

            return (poles, zeros);
        }

        private static (int N, double EpsP, double Wp, double kW) GetProperties(double fp, double fs, double dt, double Gp, double Gs)
        {
            if (!(fp < fs)) throw new InvalidOperationException("Частота пропускания должна быть меньше частоты подавления");
            if (!(fp < 1 / (2 * dt))) throw new InvalidOperationException();

            var Rp = -Gp.In_dB();
            var Rs = -Gs.In_dB();

            var eps_p = (10d.Pow(Rp / 10) - 1).Sqrt();
            var eps_s = (10d.Pow(Rs / 10) - 1).Sqrt();

            var Fp = ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
            var Fs = ToAnalogFrequency(fs, dt);  // Частота подавления аналогового пропита

            var Wp = Consts.pi2 * Fp;

            var k_eps = eps_s / eps_p;
            var k_W = Fs / Fp;

            var N = (int)Math.Ceiling(arch(k_eps) / arch(k_W)); // Порядок фильтра

            return (N, eps_p, Wp, k_W);
        }

        private static (double[] A, double[] B) InitializeI(double fp, double fs, double dt, double Gp, double Gs)
        {
            var properties = GetProperties(fp, fs, dt, Gp, Gs);
            var poles = GetAnalogPolesI(properties.N, properties.EpsP);
            var translated_poles = ToZArray(poles, dt, properties.Wp);

            return (Array.Empty<double>(), Array.Empty<double>());
        }

        private static (double[] A, double[] B) InitializeII(double fp, double fs, double dt, double Gp, double Gs)
        {
            var properties = GetProperties(fp, fs, dt, Gp, Gs);
            var polynom = GetAnalogPolesII(properties.N, properties.EpsP);
            var z_poles = ToZArray(polynom.Poles, dt, properties.Wp);
            var z_zeros = ToZArray(polynom.Zeros, dt, properties.Wp);

            //var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
            //var A = GetCoefficientsInverted(z_poles).ToRe();

            return (Array.Empty<double>(), Array.Empty<double>());
        }

        private static (double[] A, double[] B) InitializeIICorrected(double fp, double fs, double dt, double Gp, double Gs)
        {
            var properties = GetProperties(fp, fs, dt, Gp, Gs);
            var polynom = GetAnalogPolesII(properties.N, properties.EpsP);
            var translated_poles = ToZArray(polynom.Poles, dt, properties.Wp / properties.kW);
            var translated_zeros = ToZArray(polynom.Zeros, dt, properties.Wp / properties.kW);

            return (Array.Empty<double>(), Array.Empty<double>());
        }

        public ChebyshevType FilterType { get; }

        /// <summary>Инициализация нового фильтра Чебышева нижних частот</summary>
        /// <param name="fp">Частота пропускания</param>
        /// <param name="fs">Частота заграждения</param>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
        /// <param name="Gs">Затухание в полосе заграждения (0.005623413 = -45 дБ)</param>
        /// <param name="Type">Тип (род) фильтра чебышева</param>
        public ChebyshevLowPass(double fp, double fs, double dt, double Gp = 0.891250938, double Gs = 0.005623413, ChebyshevType Type = ChebyshevType.I)
            : this(Type switch
            {
                ChebyshevType.I => InitializeI(fp, fs, dt, Gp, Gs),
                ChebyshevType.II => InitializeII(fp, fs, dt, Gp, Gs),
                ChebyshevType.IICorrected => InitializeIICorrected(fp, fs, dt, Gp, Gs),
                _ => throw new InvalidEnumArgumentException(nameof(Type), (int)Type, typeof(ChebyshevType))
            }) => FilterType = Type;

        private ChebyshevLowPass((double[] A, double[] B) config) : base(config.B, config.A) { }
    }
}
