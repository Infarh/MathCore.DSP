using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;

namespace MathCore.DSP.Fourier
{
    /// <summary>Класс преобразования Фурье</summary>
    public static class FT
    {
        /* -------------------------------------------------------------------------------------------- */

        /// <summary>Коэффициент прямого преобразования Фурье</summary>
        /// <param name="N">Число спектральных составляющих</param>
        /// <returns>Комплексное значение коэффициента Фурье</returns>
        public static Complex W(int N) => Complex.Exp(1d / N, -Consts.pi2 / N);

        /// <summary>Коэффициент обратного преобразования Фурье</summary>
        /// <param name="N">Число спектральных составляющих</param>
        /// <returns>Комплексное значение обратного коэффициента Фурье</returns>
        public static Complex Winv(int N) => Complex.Exp(Consts.pi2 / N);

        /// <summary>Коэффициент прямого преобразования Фурье для спектральной составляющей</summary>
        /// <param name="k">Номер спектральной составляющей</param>
        /// <param name="N">Число спектральных составляющих</param>
        /// <returns>Комплексное значение коэффициента Фурье</returns>
        public static Complex W(int k, int N) => Complex.Exp(1d / N, -Consts.pi2 * k / N);

        /// <summary>Коэффициент обратного преобразования Фурье для спектральной составляющей</summary>
        /// <param name="k">Номер спектральной составляющей</param>
        /// <param name="N">Число спектральных составляющих</param>
        /// <returns>Комплексное значение обратного коэффициента Фурье</returns>
        public static Complex Winv(int k, int N) => Complex.Exp(1d / N, -Consts.pi2 * k / N);

        /* -------------------------------------------------------------------------------------------- */

        [NotNull]
        public static Complex[] FourierTransform([NotNull] this IEnumerable<double> Values, bool IsInverse = false)
            => Values.ToArray().FourierTransform(IsInverse);


        /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
        /// <param name="Values">Массив отсчётов функции</param>
        /// <param name="IsInverse">Обратное преобразование</param>
        [NotNull]
        public static Complex[] FourierTransform([NotNull] this double[] Values, bool IsInverse = false)
        {
            if (Values is null) throw new ArgumentNullException(nameof(Values));

            var N = Values.Length;
            var spectrum = new Complex[N];
            var w = Exp.GetCoefficients(N, IsInverse);

            for (var m = 0; m < N; m++)
            {
                var p = 0.0;
                var q = 0.0;
                for (var n = 0; n < N; n++)
                {
                    var v = Values[n];
                    var i = m * n % N;
                    var ww = w[i];
                    p += v * ww._Cos;
                    q += v * ww._Sin;
                }
                spectrum[m] = new Complex(p / N, q / N);
            }

            return spectrum;
        }

        /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
        /// <param name="Values">Массив отсчётов функции</param>
        /// <param name="IsInverse">Обратное преобразование</param>
        /// <param name="progress">Метод информирования о прогрессе операции</param>
        [NotNull]
        public static Complex[] FourierTransform([NotNull] this double[] Values, bool IsInverse, [CanBeNull] Action<double> progress)
        {
            if (Values is null) throw new ArgumentNullException(nameof(Values));

            var N = Values.Length;
            var spectrum = new Complex[N];
            var w = Exp.GetCoefficients(N, IsInverse);

            for (var m = 0; m < N; m++)
            {
                var p = 0.0;
                var q = 0.0;
                for (var n = 0; n < N; n++)
                {
                    var v = Values[n];
                    var i = m * n % N;
                    var ww = w[i];
                    p += v * ww._Cos;
                    q += v * ww._Sin;
                }
                spectrum[m] = new Complex(p / N, q / N);
                progress?.Invoke((double)m / N);
            }
            progress?.Invoke(1);

            return spectrum;
        }

        /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
        /// <param name="Values">Массив отсчётов функции</param>
        /// <param name="Inverse">Обратное преобразование</param>
        [NotNull]
        public static Complex[] FourierTransform([NotNull] this Complex[] Values, bool Inverse = false)
        {
            if (Values is null) throw new ArgumentNullException(nameof(Values));

            var spectrum = new Complex[Values.Length];
            var N = spectrum.Length;
            var w = Exp.GetCoefficients(N, Inverse);

            for (var m = 0; m < N; m++)
            {
                var P = 0.0;
                var Q = 0.0;
                for (var n = 0; n < N; n++)
                {
                    var v = Values[n];
                    var re = v.Re;
                    var im = v.Im;
                    var i = n * m % N;
                    var ww = w[i];
                    P += re * ww._Cos - im * ww._Sin;
                    Q += im * ww._Cos + re * ww._Sin;
                }
                spectrum[m] = Inverse ? new Complex(P, Q) : new Complex(P / N, Q / N);
            }
            return spectrum;
        }

        /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
        /// <param name="Values">Массив отсчётов функции</param>
        /// <param name="Inverse">Обратное преобразование</param>
        [NotNull]
        public static Complex[] FourierTransform([NotNull] this Complex[] Values, bool Inverse, [CanBeNull] Action<double> progress)
        {
            if (Values is null) throw new ArgumentNullException(nameof(Values));

            var spectrum = new Complex[Values.Length];
            var N = spectrum.Length;
            var w = Exp.GetCoefficients(N, Inverse);

            for (var m = 0; m < N; m++)
            {
                var P = 0.0;
                var Q = 0.0;
                for (var n = 0; n < N; n++)
                {
                    var v = Values[n];
                    var re = v.Re;
                    var im = v.Im;
                    var i = n * m % N;
                    var ww = w[i];
                    P += re * ww._Cos - im * ww._Sin;
                    Q += im * ww._Cos + re * ww._Sin;
                }

                progress?.Invoke((double)m / N);
                spectrum[m] = Inverse ? new Complex(P, Q) : new Complex(P / N, Q / N);
            }
            progress?.Invoke(1);
            return spectrum;
        }
    }
}
