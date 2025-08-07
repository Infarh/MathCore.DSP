using static System.Math;

using static MathCore.Consts;

namespace MathCore.DSP.Fourier;

/// <summary>Быстрое преобразование Фурье</summary>
public static class FFT
{
    /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
    /// <param name="Values">Массив отсчётов функции</param>
    public static Complex[] FastFourierTransform(this double[] Values)
    {
        return Fourier.fft.FFT(Values);

        //var N = Values.Length;
        //if(!N.IsPowerOf2()) N = 1 << ((int)Log(N, 2) + 1);

        //var real_spectrum = new double[N * 2];
        //N = Values.Length;
        //for(var n = 0; n < N; n++)
        //    real_spectrum[2 * n] = Values[n];

        //fft(ref real_spectrum, false);

        //var spectrum = new Complex[Values.Length];
        //for(int m = 0, M = spectrum.Length; m < M; m++)
        //    spectrum[m] = new Complex(real_spectrum[2 * m], real_spectrum[2 * m + 1]);

        //return spectrum;
    }

    /// <summary>Прямое преобразование отсчётов функции в спектр</summary>
    /// <param name="Values">Массив отсчётов функции</param>
    public static Complex[] FastFourierTransform(this Complex[] Values)
    {
        return Fourier.fft.FFT(Values);

        //var length = Values.Length;
        //if(!length.IsPowerOf2())
        //{
        //    var length_log2 = Log(length, 2);
        //    length_log2 -= Round(length_log2);
        //    length += (int)Pow(2, length_log2);
        //}

        //var spectrum = new double[length * 2];
        //var values_length = Values.Length;
        //for(var i = 0; i < values_length; i++)
        //{
        //    spectrum[2 * i] = Values[i].Re;
        //    spectrum[2 * i + 1] = Values[i].Im;
        //}

        //fft(ref spectrum, false);

        //var Spectrum = new Complex[length];
        //values_length = Spectrum.Length;
        //for(var i = 0; i < values_length; i++)
        //    Spectrum[i] = new Complex(spectrum[2 * i], spectrum[2 * i + 1]);


        //return Spectrum;
    }

    /// <summary>Обратное преобразование отсчётов спектра в отсчёты сигнала</summary>
    /// <param name="Spectrum">Массив отсчётов спектра</param>
    public static Complex[] FastFourierInverse(this Complex[] Spectrum)
    {
        return Fourier.fft.FFT_Complex_Inverse(Spectrum);

        //var spectrum_length = Spectrum.Length;
        //if (!spectrum_length.IsPowerOf2())
        //{

        //    var spectrum_length_log2 = Log(spectrum_length, 2);
        //    spectrum_length = 1 << (1 + (int)Math.Floor(spectrum_length_log2));
        //}

        //var spectrum = new double[spectrum_length * 2];
        //var values_length = Spectrum.Length;
        //for (var i = 0; i < values_length; i++)
        //{
        //    spectrum[2 * i] = Spectrum[i].Re;
        //    spectrum[2 * i + 1] = Spectrum[i].Im;
        //}

        //fft(ref spectrum, false);

        //var samples = new Complex[spectrum_length];
        //values_length = samples.Length;
        //for (var i = 0; i < values_length; i++)
        //    samples[i] = new Complex(spectrum[2 * i], spectrum[2 * i + 1]);

        //return samples.ResamplingOptimal(Spectrum.Length);
    }

    // ReSharper disable once InconsistentNaming
    private static void fft(ref double[] Values, bool IsInverse)
    {
        var N = Values.Length / 2;
        if (!N.IsPowerOf2())
            throw new ArgumentException("Число элементов выборки должно быть степенью двойки");

        var i_sign = IsInverse ? -1 : 1;

        var n = N << 1;
        var j = 1;

        //Операция бабочка
        for (var ii = 1; ii <= N; ii++)
        {
            var i = (ii << 1) - 1;
            if (j > i)
            {
                var temp_r = Values[j - 1];
                var temp_i = Values[j];
                Values[j - 1] = Values[i - 1];
                Values[j] = Values[i];
                Values[i - 1] = temp_r;
                Values[i] = temp_i;
            }
            var m = n >> 1;
            while (m >= 2 && j > m)
            {
                j -= m;
                m >>= 1;
            }
            j += m;
        }

        var m_max = 2;
        var theta0 = pi2 * i_sign;
        while (n > m_max)
        {
            var i_step = m_max << 1;
            var theta = theta0 / m_max;
            var w_pr = Sin(.5 * theta);
            w_pr *= -2 * w_pr;
            var w_pi = Sin(theta);
            double w_r = 1;
            double w_i = 0;
            for (var ii = 1; ii <= (m_max >> 1); ii++)
            {
                var m = (ii << 1) - 1;
                for (var jj = 0; jj <= (n - m) / i_step; jj++)
                {
                    var i = m + jj * i_step;
                    j = i + m_max;
                    var temp_r = w_r * Values[j - 1] - w_i * Values[j];
                    var temp_i = w_r * Values[j] + w_i * Values[j - 1];
                    Values[j - 1] = Values[i - 1] - temp_r;
                    Values[j] = Values[i] - temp_i;
                    Values[i - 1] += temp_r;
                    Values[i] += temp_i;
                }
                var w_temp = w_r;
                w_r = w_r * w_pr - w_i * w_pi + w_r;
                w_i = w_i * w_pr + w_temp * w_pi + w_i;
            }
            m_max = i_step;
        }

        if (IsInverse) return;
        for (var i = 1; i <= 2 * N; i++)
            Values[i - 1] /= N;
    }

    /// <summary>Целочисленное преобразование Фурье</summary>
    /// <param name="a">Массив целых чисел</param>
    /// <param name="invert">Обратное преобразование</param>
    public static void FFT_int(this int[] a, bool invert = false)
    {
        var n = a.Length;

        for (int i = 1, j = 0; i < n; i++)
        {
            var bit = n >> 1;
            for (; j >= bit; bit >>= 1)
                j -= bit;

            j += bit;

            if (i >= j) continue;
            a[i] ^= a[j];
            a[j] ^= a[i];
            a[i] ^= a[j];
        }

        const int root = 5;
        const int root_1 = 4404020;
        for (var len = 2; len <= n; len <<= 1)
        {
            var w_len = invert ? root_1 : root;

            const int root_pw = 1 << 20;
            const int mod = 7340033;

            for (var i = len; i < root_pw; i <<= 1)
                w_len = w_len * 1 * w_len % mod;
            for (var i = 0; i < n; i += len)
            {
                var w = 1;
                for (var j = 0; j < len / 2; ++j)
                {
                    var u = a[i + j];
                    var v = a[i + j + len >> 1] * 1 * w % mod;
                    a[i + j] = u + v < mod ? u + v : u + v - mod;
                    a[i + j + len / 2] = u - v >= 0 ? u - v : u - v + mod;
                    w = w * 1 * w_len % mod;
                }
            }

            if (!invert) continue;
            var nrev = ReverseMod(n, mod);
            for (var i = 0; i < n; i++)
                a[i] = a[i] * 1 * nrev % mod;
        }
    }

    private static int EuclidEx(int n, int mod, out int x, out int y)
    {
        if (mod == 0)
        {
            x = 1;
            y = 0;
            return n;
        }
        var x2 = 1;
        var x1 = 0;
        var y2 = 0;
        var y1 = 1;
        while (mod > 0)
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

    public static double[] Recursive_FFT(double[] a)
    {
        var n = a.Length;
        return n.IsPowerOf2() ? Recursive_FFTInternal(a, n) : throw new ArgumentException("Длина массива должна быть степенью 2", nameof(a));
    }

    private static double[] Recursive_FFTInternal(double[] a, int n)
    {
        switch (n)
        {
            case 1: return a;
            case 2: return [a[0] + a[1], a[0] - a[1]];
        }

        var n2 = n / 2;

        var even = new double[n2];
        var odd = new double[n2];
        for (var i = 0; i < n / 2; i++)
        {
            even[i] = a[i * 2];
            odd[i] = a[i * 2 + 1];
        }

        even = Recursive_FFTInternal(even, n2);
        odd = Recursive_FFTInternal(odd, n2);

        var y = new double[n];
        var wn = Exp(pi2 / n);
        var w = 1d;
        for (var k = 0; k < n2; k++)
        {
            y[k] = even[k] + w * odd[k];
            y[k + n2] = even[k] - w * odd[k];
            w *= wn;
        }
        return y;
    }

    public static Task<double[]> Recursive_FFTAsync(double[] a, int MinAsyncLength = 256)
    {
        var n = a.Length;
        return n.IsPowerOf2()
            ? Recursive_FFTInternalAsync(a, n, MinAsyncLength)
            : throw new ArgumentException("Длина массива должна быть степенью 2", nameof(a));
    }

    private static async Task<double[]> Recursive_FFTInternalAsync(double[] a, int n, int MinAsyncLength = 256)
    {
        switch (n)
        {
            case 1:
                return a;
            case 2:
                return [a[0] + a[1], a[0] - a[1]];
        }
        var n2 = n / 2;

        var even = new double[n2];
        var odd = new double[n2];
        for (var i = 0; i < n / 2; i++)
        {
            even[i] = a[i * 2];
            odd[i] = a[i * 2 + 1];
        }

        if (n <= MinAsyncLength)
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
        var wn = Exp(pi2 / n);
        var w = 1d;
        for (var k = 0; k < n2; k++)
        {
            y[k] = even[k] + w * odd[k];
            y[k + n2] = even[k] - w * odd[k];
            w *= wn;
        }
        return y;
    }

    /// <summary>Быстрое преобразование Фурье по алгоритму Блюштейна (Bluestein) для произвольной длины</summary>
    /// <param name="Values">Массив комплексных отсчётов</param>
    /// <returns>Массив комплексных коэффициентов спектра</returns>
    public static Complex[] FFT_Bluestein(this Complex[] Values)
    {
        // Алгоритм Блюштейна позволяет вычислять БПФ для любого N через свёртку
        var n = Values.Length;
        if (n == 0) return [];
        if (n == 1) return [Values[0]];
        if (n.IsPowerOf2()) return Values.FastFourierTransform();

        // Находим минимальную степень двойки >= 2n-1
        var m = 1;
        while (m < 2 * n - 1) m <<= 1;

        var a = new Complex[m]; // входной вектор
        var b = new Complex[m]; // вектор для свёртки
        //var w = Complex.Exp(-pi2 / n); // корень из единицы

        // Предвычисляем фазовые множители
        for (var k = 0; k < n; k++)
        {
            var angle = (k * k) % (2 * n);
            var wk = Complex.Exp(pi2 * angle / n);
            a[k] = Values[k] * wk; // a[k] = x[k] * w^{k^2}
            b[k] = Complex.Exp(-pi2 * angle / n); // b[k] = w^{-k^2}
        }

        Array.Clear(a, n, m - n); //for (var k = n; k < m; k++) b[k] = Complex.Zero;
        Array.Clear(b, n, m - n); //for (var k = n; k < m; k++) a[k] = Complex.Zero;

        // b[m - k] = b[k] для k=1..n-1 (симметрия)
        for (var k = 1; k < n; k++) b[m - k] = b[k];

        // Выполняем свёртку через БПФ
        var A = a.FastFourierTransform();
        var B = b.FastFourierTransform();
        var C = new Complex[m];
        for (var i = 0; i < m; i++) C[i] = A[i] * B[i];
        var c = C.FastFourierInverse();

        // Окончательный результат
        var result = new Complex[n];
        for (var k = 0; k < n; k++)
        {
            var angle = (k * k) % (2 * n);
            var wk = Complex.Exp(pi2 * angle / n);
            result[k] = c[k] * wk;
        }
        return result;
    }

    /// <summary>Быстрое преобразование Фурье по алгоритму Радера (Rader) для простых длин</summary>
    /// <param name="Values">Массив комплексных отсчётов</param>
    /// <returns>Массив комплексных коэффициентов спектра</returns>
    public static Complex[] FFT_Rader(this Complex[] Values)
    {
        // Алгоритм Радера работает для простых длин
        var n = Values.Length;
        if (n == 0) return [];
        if (n == 1) return [Values[0]];
        if (n.IsPowerOf2()) return Values.FastFourierTransform();
        if (!IsPrime(n)) throw new ArgumentException("Длина массива должна быть простым числом", nameof(Values));

        // Находим примитивный корень g по модулю n
        var g = FindPrimitiveRoot(n);
        var g_inv = ModInverse(g, n);

        var a = new Complex[n - 1];
        var b = new Complex[n - 1];
        for (var i = 0; i < n - 1; i++)
        {
            var j = ModPow(g, i, n);
            a[i] = Values[j];
            b[i] = Complex.Exp(-pi2 * j / n);
        }

        // Дополняем до степени двойки
        var m = 1;
        while (m < 2 * (n - 1)) m <<= 1;
        var a_pad = new Complex[m];
        var b_pad = new Complex[m];
        for (var i = 0; i < n - 1; i++)
        {
            a_pad[i] = a[i];
            b_pad[i] = b[i];
        }
        for (var i = n - 1; i < m; i++)
        {
            a_pad[i] = Complex.Zero;
            b_pad[i] = Complex.Zero;
        }

        // Свёртка через БПФ
        var A = a_pad.FastFourierTransform();
        var B = b_pad.FastFourierTransform();
        var C = new Complex[m];
        for (var i = 0; i < m; i++) C[i] = A[i] * B[i];
        var c = C.FastFourierInverse();

        // Формируем результат
        var result = new Complex[n];
        var sum = Values[0];
        for (var i = 0; i < n; i++) sum += Values[i];
        result[0] = sum;
        for (var k = 1; k < n; k++)
        {
            var idx = ModPow(g_inv, k - 1, n);
            result[k] = Values[0] + c[idx];
        }
        return result;
    }

    /// <summary>Проверка простоты числа</summary>
    private static bool IsPrime(int n)
    {
        if (n < 2) return false;
        if (n == 2) return true;
        if (n % 2 == 0) return false;
        for (var i = 3; i * i <= n; i += 2)
            if (n % i == 0) return false;
        return true;
    }

    /// <summary>Возведение в степень по модулю</summary>
    private static int ModPow(int a, int exp, int mod)
    {
        var res = 1;
        for (; exp > 0; exp >>= 1, a = a * a % mod)
            if ((exp & 1) != 0) res = res * a % mod;
        return res;
    }

    /// <summary>Обратный элемент по модулю</summary>
    private static int ModInverse(int a, int mod)
    {
        var t = 0; var newt = 1;
        var r = mod; var newr = a;
        while (newr != 0)
        {
            var quotient = r / newr;
            (t, newt) = (newt, t - quotient * newt);
            (r, newr) = (newr, r - quotient * newr);
        }
        if (r > 1) throw new ArgumentException("Число не обратимо");
        if (t < 0) t += mod;
        return t;
    }

    /// <summary>Поиск примитивного корня по модулю простого числа</summary>
    private static int FindPrimitiveRoot(int p)
    {
        var phi = p - 1;
        var factors = new List<int>();
        var n = phi;
        for (var i = 2; i * i <= n; i++)
            if (n % i == 0)
            {
                factors.Add(i);
                while (n % i == 0) n /= i;
            }
        if (n > 1) factors.Add(n);
        for (var res = 2; res <= p; res++)
        {
            var ok = true;
            foreach (var factor in factors)
                if (ModPow(res, phi / factor, p) == 1)
                {
                    ok = false;
                    break;
                }
            if (ok) return res;
        }
        throw new ArgumentException("Не найден примитивный корень");
    }
}