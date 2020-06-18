using System;
using System.Threading.Tasks;
using MathCore.Annotations;

namespace MathCore.DSP.Fourier
{
    /// <summary>Быстрое преобразование Фурье</summary>
    public static class FFT
    {
        /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
        /// <param name="Values">Массив отсчётов функции</param>
        [NotNull]
        public static Complex[] FastFourierTransform([NotNull] this double[] Values)
        {
            var N = Values.Length;
            if(!N.IsPowerOf2()) N = 1 << ((int)Math.Log(N, 2) + 1);

            var real_spectrum = new double[N * 2];
            N = Values.Length;
            for(var n = 0; n < N; n++)
                real_spectrum[2 * n] = Values[n];

            fft(ref real_spectrum, false);

            var spectrum = new Complex[Values.Length];
            for(int m = 0, M = spectrum.Length; m < M; m++)
                spectrum[m] = new Complex(real_spectrum[2 * m], real_spectrum[2 * m + 1]);

            return spectrum;
        }

        /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
        /// <param name="Values">Массив отсчётов функции</param>
        [NotNull]
        public static Complex[] FastFourierTransform([NotNull] this Complex[] Values)
        {
            var length = Values.Length;
            if(!length.IsPowerOf2())
            {
                var length_log2 = Math.Log(length, 2);
                length_log2 -= Math.Round(length_log2);
                length += (int)Math.Pow(2, length_log2);
            }

            var spectrum = new double[length * 2];
            var values_length = Values.Length;
            for(var i = 0; i < values_length; i++)
            {
                spectrum[2 * i] = Values[i].Re;
                spectrum[2 * i + 1] = Values[i].Im;
            }

            fft(ref spectrum, false);

            var Spectrum = new Complex[length];
            values_length = Spectrum.Length;
            for(var i = 0; i < values_length; i++)
                Spectrum[i] = new Complex(spectrum[2 * i], spectrum[2 * i + 1]);


            return Spectrum;
        }

        /// <summary>Обратное преобразование отсчётов спектра в отсчёты сигнала</summary>
        /// <param name="Spectrum">Массив отсчётов спектра</param>
        [NotNull]
        public static Complex[] FastFurierInverse([NotNull] this Complex[] Spectrum)
        {
            var spectrum_length = Spectrum.Length;
            if(!spectrum_length.IsPowerOf2())
            {
                var spectrum_length_log2 = Math.Log(spectrum_length, 2);
                spectrum_length_log2 -= Math.Round(spectrum_length_log2);
                spectrum_length += (int)Math.Pow(2, spectrum_length_log2);
            }

            var values = new double[spectrum_length * 2];
            for(var i = 0; i < spectrum_length; i++)
            {
                values[2 * i] = Spectrum[i].Re;
                values[2 * i + 1] = Spectrum[i].Im;
            }

            fft(ref values, true);

            var complex_values = new Complex[spectrum_length];
            var complex_values_length = complex_values.Length;
            for(var i = 0; i < complex_values_length; i++)
                complex_values[i] = new Complex(values[2 * i], values[2 * i + 1]);

            return complex_values;
        }

#pragma warning disable IDE1006 // Стили именования
        private static void fft([NotNull] ref double[] Values, bool IsInverse)
#pragma warning restore IDE1006 // Стили именования
        {
            var N = Values.Length / 2;
            if(!N.IsPowerOf2())
                throw new ArgumentException("Число элементов выборки должно быть степенью двойки");

            var i_sign = IsInverse ? -1 : 1;

            var n = N << 1;
            var j = 1;

            //Операция бабочка
            for(var ii = 1; ii <= N; ii++)
            {
                var i = (ii << 1) - 1;
                if(j > i)
                {
                    var temp_r = Values[j - 1];
                    var temp_i = Values[j];
                    Values[j - 1] = Values[i - 1];
                    Values[j] = Values[i];
                    Values[i - 1] = temp_r;
                    Values[i] = temp_i;
                }
                var m = n >> 1;
                while(m >= 2 && j > m)
                {
                    j -= m;
                    m >>= 1;
                }
                j += m;
            }

            var m_max = 2;
            var theta0 = Consts.pi2 * i_sign;
            while(n > m_max)
            {
                var i_step = m_max << 1;
                var theta = theta0 / m_max;
                var w_pr = Math.Sin(.5 * theta);
                w_pr *= -2 * w_pr;
                var w_pi = Math.Sin(theta);
                double w_r = 1;
                double w_i = 0;
                for(var ii = 1; ii <= (m_max >> 1); ii++)
                {
                    var m = (ii << 1) - 1;
                    for(var jj = 0; jj <= (n - m) / i_step; jj++)
                    {
                        var i = m + jj * i_step;
                        j = i + m_max;
                        var temp_r = w_r * Values[j - 1] - w_i * Values[j];
                        var temp_i = w_r * Values[j] + w_i * Values[j - 1];
                        Values[j - 1] = Values[i - 1] - temp_r;
                        Values[j] = Values[i] - temp_i;
                        Values[i - 1] = Values[i - 1] + temp_r;
                        Values[i] = Values[i] + temp_i;
                    }
                    var w_temp = w_r;
                    w_r = w_r * w_pr - w_i * w_pi + w_r;
                    w_i = w_i * w_pr + w_temp * w_pi + w_i;
                }
                m_max = i_step;
            }

            if(IsInverse) return;
            for(var i = 1; i <= 2 * N; i++)
                Values[i - 1] = Values[i - 1] / N;
        }

        /// <summary>Целочисленное преобразование Фурье</summary>
        /// <param name="a">Массив целых чисел</param>
        /// <param name="invert">Обратное преобразование</param>
        public static void FFT_int([NotNull] this int[] a, bool invert = false)
        {
            var n = a.Length;

            for(int i = 1, j = 0; i < n; i++)
            {
                var bit = n >> 1;
                for(; j >= bit; bit >>= 1)
                    j -= bit;

                j += bit;

                if(i >= j) continue;
                a[i] ^= a[j];
                a[j] ^= a[i];
                a[i] ^= a[j];
            }

            const int root = 5;
            const int root_1 = 4404020;
            for(var len = 2; len <= n; len <<= 1)
            {
                var wlen = invert ? root_1 : root;

                const int root_pw = 1 << 20;
                const int mod = 7340033;

                for(var i = len; i < root_pw; i <<= 1)
                    wlen = wlen * 1 * wlen % mod;
                for(var i = 0; i < n; i += len)
                {
                    var w = 1;
                    for(var j = 0; j < len / 2; ++j)
                    {
                        var u = a[i + j];
                        var v = a[i + j + len >> 1] * 1 * w % mod;
                        a[i + j] = u + v < mod ? u + v : u + v - mod;
                        a[i + j + len / 2] = u - v >= 0 ? u - v : u - v + mod;
                        w = w * 1 * wlen % mod;
                    }
                }

                if(!invert) continue;
                var nrev = ReverseMod(n, mod);
                for(var i = 0; i < n; i++)
                    a[i] = a[i] * 1 * nrev % mod;
            }
        }

        private static int EuclidEx(int n, int mod, out int x, out int y)
        {
            if(mod == 0)
            {
                x = 1;
                y = 0;
                return n;
            }
            var x2 = 1;
            var x1 = 0;
            var y2 = 0;
            var y1 = 1;
            while(mod > 0)
            {
                var q = n / mod;
                var r = n - q * mod;
                x = x2 - q * x1;
                y = y2 - q * y1;
                n = mod;
                mod = r;
                x2 = x1;
                x1 = x;
                y2 = y1;
                y1 = y;
            }
            x = x2;
            y = y2;
            return n;
        }

        private static int ReverseMod(int n, int mod)
        {
            var d = EuclidEx(n, mod, out var x, out _);
            return d == 1 ? x : 0;
        }


        public static double[] Recursive_FFT([NotNull] double[] a)
        {
            var n = a.Length;
            if(!n.IsPowerOf2())
                throw new ArgumentException("Длина массива должна быть степенью 2", nameof(a));
            return Recursive_FFTInternal(a, n);
        }
        private static double[] Recursive_FFTInternal(double[] a, int n)
        {
            switch (n)
            {
                case 1: return a;
                case 2: return new[] { a[0] + a[1], a[0] - a[1] };
            }

            var n2 = n / 2;

            var even = new double[n2];
            var odd = new double[n2];
            for(var i = 0; i < n / 2; i++)
            {
                even[i] = a[i * 2];
                odd[i] = a[i * 2 + 1];
            }

            even = Recursive_FFTInternal(even, n2);
            odd = Recursive_FFTInternal(odd, n2);

            var y = new double[n];
            var wn = Math.Exp(Consts.pi2 / n);
            var w = 1d;
            for(var k = 0; k < n2; k++)
            {
                y[k] = even[k] + w * odd[k];
                y[k + n2] = even[k] - w * odd[k];
                w *= wn;
            }
            return y;
        }

        [NotNull, ItemNotNull]
        public static Task<double[]> Recursive_FFTAsync([NotNull] double[] a, int MinAsyncLength = 256)
        {
            var n = a.Length;
            if(!n.IsPowerOf2())
                throw new ArgumentException("Длина массива должна быть степенью 2", nameof(a));
            return Recursive_FFTInternalAsync(a, n, MinAsyncLength);
        }

        [ItemNotNull]
        private static async Task<double[]> Recursive_FFTInternalAsync([NotNull] double[] a, int n, int MinAsyncLength = 256)
        {
            switch(n)
            {
                case 1:
                    return a;
                case 2:
                    return new[] { a[0] + a[1], a[0] - a[1] };
            }
            var n2 = n / 2;

            var even = new double[n2];
            var odd = new double[n2];
            for(var i = 0; i < n / 2; i++)
            {
                even[i] = a[i * 2];
                odd[i] = a[i * 2 + 1];
            }

            if(n <= MinAsyncLength)
            {
                even = Recursive_FFTInternal(even, n2);
                odd = Recursive_FFTInternal(odd, n2);
            }
            else
            {
                var t_even = Recursive_FFTInternalAsync(even, n2, MinAsyncLength);
                var t_odd = Recursive_FFTInternalAsync(odd, n2, MinAsyncLength);
                await Task.WhenAll(t_even, t_odd).ConfigureAwait(false);
                even = t_even.Result;
                odd = t_odd.Result;
            }

            var y = new double[n];
            var wn = Math.Exp(Consts.pi2 / n);
            var w = 1d;
            for(var k = 0; k < n2; k++)
            {
                y[k] = even[k] + w * odd[k];
                y[k + n2] = even[k] - w * odd[k];
                w *= wn;
            }
            return y;
        }


        //[Copyright("", url = "http://rain.ifmo.ru/cat/view.php/theory/math/fft-2004")]
        //public static double[] Recursive_FFT(double[] a)
        //{
        //    var n = a.Length;
        //    if(n == 1) return a;
        //    var wn = Complex.Exp(Consts.pi2/n);
        //    var w = 1;
        //    var y1 = Recursive_FFT()
        //}
    }
}
