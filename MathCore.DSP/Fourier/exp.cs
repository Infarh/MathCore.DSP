
using System;
using MathCore.Annotations;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Fourier
{
    /// <summary>Комплексное экспоненциальное значение</summary>
    internal struct Exp
    {
        /// <summary>Синусная составляющая</summary>
        public readonly double Sin;

        /// <summary>Косинусная составляющая</summary>
        public readonly double Cos;

        /// <summary>Инициализация нового комплексного экспоненциального значения</summary>
        /// <param name="Sin">Синусная составляющая</param>
        /// <param name="Cos">Косинусная составляющая</param>
        public Exp(double Sin, double Cos) { this.Sin = Sin; this.Cos = Cos; }

        /// <summary>Рассчитать массив комплексных экспоненциальных коэффициентов</summary>
        /// <param name="N"></param>
        /// <param name="IsInverse"></param>
        /// <returns></returns>
        [NotNull]
        public static Exp[] GetCoefficients(int N, bool IsInverse = false)
        {
            var w = new Exp[N];
            var darg = Consts.pi2 / N;
            if(IsInverse) darg *= -1;
            var arg = 0.0;
            for(var i = 0; i < N; i++)
            {
                w[i] = new Exp(Math.Sin(arg), Math.Cos(arg));
                arg += darg;
            }
            return w;
        }
    }
}