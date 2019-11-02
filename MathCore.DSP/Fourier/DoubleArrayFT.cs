using System;
using MathCore.Annotations;

namespace MathCore.DSP.Fourier
{
    using Spectrum = Func<int, Complex>;

    /// <summary>Методы-расширения для вещественного и комплексного массивов</summary>
    public static class DoubleArrayFT
    {
        /// <summary>Выполнить преобразование Фурье для вещественного массива</summary>
        /// <param name="s">Массив вещественных значений отсчётов</param>
        /// <param name="IsInverse">Выполнить обратное преобразование</param>
        /// <returns>Спектр</returns>
        [NotNull]
        public static Spectrum GetFourierTransformation([NotNull] this double[] s, bool IsInverse = false)
        {
            var N = s.Length;
            var w = Exp.GetCoefficients(N, IsInverse);

            return m =>
            {
                var P = 0.0;
                var Q = 0.0;
                for(var n = 0; n < N; n++)
                {
                    var val = s[n];
                    var ww = w[n * m % N];
                    P += val * ww._Cos;
                    Q += val * ww._Sin;
                }

                return new Complex(P / N, Q / N);
            };
        }

        /// <summary>Выполнить преобразование Фурье для комплексного массива</summary>
        /// <param name="s">Массив комплексных значений отсчётов</param>
        /// <param name="IsInverse">Выполнить обратное преобразование</param>
        /// <returns>Спектр</returns>
        [NotNull]
        public static Spectrum GetFourierTransformation([NotNull] this Complex[] s, bool IsInverse = false)
        {
            var N = s.Length;
            var w = Exp.GetCoefficients(N, IsInverse);

            return m =>
            {
                var P = 0.0;
                var Q = 0.0;
                for(var n = 0; n < N; n++)
                {
                    var (re, im) = s[n];
                    var ww = w[n * m % N];
                    var sin = ww._Sin;
                    var cos = ww._Cos;
                    P += re * cos - im * sin;
                    Q += im * cos + re * sin;
                }

                return new Complex(P / N, Q / N);
            };
        }
    }
}