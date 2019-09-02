﻿using System;

namespace MathCore.DSP.Fourier
{
    using IntSpectrum = Func<int, Complex>;
    using DoubleSpectrum = Func<double, Complex>;

    public static class DoubleFuncFT
    {
        /// <summary>Выполнить преобразвоание Фурье</summary>
        /// <param name="s">Вечественная функция</param>
        /// <param name="t1">Начало интервала</param>
        /// <param name="t2">Конец интервала</param>
        /// <param name="IsInverse">Обратное преобразование</param>
        /// <param name="dt">Шаг численного расчёта</param>
        /// <returns>Спектр</returns>
        public static DoubleSpectrum GetFurierTransformation(
            this Func<double, double> s,
            double t1,
            double t2,
            bool IsInverse = false,
            double dt = 1e-4)
        {
            var delta_t = t2 - t1;
            var N = (int)Math.Abs(delta_t / dt);
            dt = delta_t / N;
            var w = Consts.pi2;
            if(IsInverse) w *= -1;

            return f =>
            {
                var pif = w * f;
                if(IsInverse) pif *= -1;
                var t = t1;
                var val = s(t);
                var arg = pif * t;
                var p = val * Math.Cos(arg);
                var q = val * Math.Sin(arg);
            
                var re_s = .0;
                var im_s = .0;
                for(var i = 0; i < N; i++)
                {
                    val = s(t);
                    arg = pif * t;
                    re_s += p + (p = val * Math.Cos(arg));
                    im_s += q + (q = val * Math.Sin(arg));
                    t += dt;
                }
                return new Complex(re_s * .5 * dt, im_s * .5 * dt);
            };
        }

        /// <summary>Выполнить преобразвоание Фурье</summary>
        /// <param name="s">Комплексная функция</param>
        /// <param name="t1">Начало интервала</param>
        /// <param name="t2">Конец интервала</param>
        /// <param name="IsInverse">Обратное преобразование</param>
        /// <param name="dt">Шаг численного расчёта</param>
        /// <returns>Спектр</returns>
        public static DoubleSpectrum GetFurierTransformation(
            this Func<double, Complex> s,
            double t1, 
            double t2,
            bool IsInverse = false,
            double dt = 1e-4)
        {
            var delta_t = t2 - t1;
            var N = (int)Math.Abs(delta_t / dt);
            dt = delta_t / N;
            var w = Consts.pi2;
            if(IsInverse) w *= -1;

            return f =>
            {   //todo: проверить алгоритм!
                var pif = w * f;
                var t = t1;
                var z = s(t).Rotate(pif * t);

                var (re, im) = z;
                for(var i = 0; i < N; i++)
                {
                    z = s(t).Rotate(pif * t);
                    re += z.Re;
                    im += z.Im;
                    t += dt;
                }
                return .5 * dt * new Complex(re, im);
            };
        }

        public static IntSpectrum GetFurierSpectrum(
            this Func<double, double> s,
            double t1,
            double t2,
            bool IsInverse = false,
            double dt = 1e-4)
        {
            var delta_t = t2 - t1;
            var N = (int)Math.Abs(delta_t / dt);
            var ss = new double[N];
            dt = delta_t / N;
            for(var n = 0; n < N; n++)
                ss[n] = s(t1 + n * dt);
            return ss.GetFourierTransformation(IsInverse);
        }

        public static IntSpectrum GetFurierSpectrum(
            this Func<double, Complex> s,
            double t1, double t2,
            bool IsInverse = false,
            double dt = 1e-4)
        {
            var delta_t = t2 - t1;
            var N = (int)Math.Abs(delta_t / dt);

            var ss = new Complex[N];
            dt = delta_t / N;
            for(var n = 0; n < N; n++)
                ss[n] = s(t1 + n * dt);
            return ss.GetFourierTransformation(IsInverse);
        }
    }
}