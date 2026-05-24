#nullable enable

using ComplexN = System.Numerics.Complex;

namespace MathCore.DSP.SDR;

/// <summary>Инструменты отладки спектра для пошагового анализа SDR-конвейера</summary>
public static class SdrDebugSpectrum
{
    /// <summary>Сформировать спектральный снимок комплексного сигнала для анализа в отладчике</summary>
    /// <param name="Samples">Массив комплексных отсчётов</param>
    /// <param name="SampleRateHz">Частота дискретизации сигнала</param>
    /// <param name="FftSize">Размер БПФ и размер выходных массивов</param>
    /// <param name="ApplyHannWindow">Признак применения окна Ханна перед БПФ</param>
    /// <returns>Кортеж из массива частот и массива уровней мощности в dB</returns>
    /// <exception cref="ArgumentNullException">Входной массив не задан</exception>
    /// <exception cref="ArgumentOutOfRangeException">Некорректно задана частота дискретизации или размер БПФ</exception>
    /// <exception cref="ArgumentException">Размер БПФ не является степенью двойки</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// #if DEBUG
    /// var (freq_hz, power_db) = SdrDebugSpectrum.GetSpectrumSnapshot(iq_data, SampleRateHz: 10_000_000, FftSize: 32768);
    /// var peak = SdrDebugSpectrum.FindPeakInRange(freq_hz, power_db, 600_000, 1_000_000);
    /// Console.WriteLine($"FFT peak: f={peak.frequency_hz:F0} Hz, p={peak.power_db:F2} dB");
    /// #endif
    /// ]]></code>
    /// </remarks>
    public static (double[] frequency_hz, double[] power_db) GetSpectrumSnapshot(
        ComplexN[] Samples,
        double SampleRateHz,
        int FftSize = 8192,
        bool ApplyHannWindow = true)
    {
        ArgumentNullException.ThrowIfNull(Samples);

        if (SampleRateHz <= 0)
            throw new ArgumentOutOfRangeException(nameof(SampleRateHz), SampleRateHz, "Частота дискретизации должна быть положительной");

        if (FftSize <= 0)
            throw new ArgumentOutOfRangeException(nameof(FftSize), FftSize, "Размер БПФ должен быть положительным");

        if ((FftSize & (FftSize - 1)) != 0)
            throw new ArgumentException("Размер БПФ должен быть степенью двойки", nameof(FftSize));

        var buffer = new ComplexN[FftSize];
        var copy_count = Math.Min(FftSize, Samples.Length);

        for (var index = 0; index < copy_count; index++)
        {
            var sample = Samples[index];
            if (!ApplyHannWindow)
            {
                buffer[index] = sample;
                continue;
            }

            var window = 0.5 - 0.5 * Math.Cos(2 * Math.PI * index / (FftSize - 1));
            buffer[index] = sample * window;
        }

        var spectrum = FastFft(buffer);
        var frequency = new double[FftSize];
        var power_db = new double[FftSize];

        for (var index = 0; index < FftSize; index++)
        {
            var shifted = (index + FftSize / 2) % FftSize;
            var frequency_hz = (index - FftSize / 2) * SampleRateHz / FftSize;
            var magnitude = spectrum[shifted].Magnitude;
            var power = magnitude * magnitude;

            frequency[index] = frequency_hz;
            power_db[index] = 10 * Math.Log10(Math.Max(power, 1e-18));
        }

        return (frequency, power_db);
    }

    /// <summary>Найти пик спектра в заданном частотном диапазоне</summary>
    /// <param name="FrequencyHz">Массив частотных бинов</param>
    /// <param name="PowerDb">Массив уровней мощности</param>
    /// <param name="FromHz">Левая граница диапазона поиска</param>
    /// <param name="ToHz">Правая граница диапазона поиска</param>
    /// <returns>Кортеж частоты и уровня найденного пика</returns>
    /// <exception cref="ArgumentNullException">Один из входных массивов не задан</exception>
    /// <exception cref="ArgumentException">Размеры массивов не совпадают</exception>
    /// <remarks>
    /// Пример использования:
    /// <code><![CDATA[
    /// var peak_zero = SdrDebugSpectrum.FindPeakInRange(freq_hz, power_db, -100_000, 100_000);
    /// var peak_rf = SdrDebugSpectrum.FindPeakInRange(freq_hz, power_db, 600_000, 1_000_000);
    /// ]]></code>
    /// </remarks>
    public static (double frequency_hz, double power_db) FindPeakInRange(double[] FrequencyHz, double[] PowerDb, double FromHz, double ToHz)
    {
        ArgumentNullException.ThrowIfNull(FrequencyHz);
        ArgumentNullException.ThrowIfNull(PowerDb);

        if (FrequencyHz.Length != PowerDb.Length)
            throw new ArgumentException("Массивы частот и мощности должны иметь одинаковую длину");

        if (FrequencyHz.Length == 0)
            return (0, double.NegativeInfinity);

        var left = Math.Min(FromHz, ToHz);
        var right = Math.Max(FromHz, ToHz);

        var best_index = -1;
        var best_power = double.NegativeInfinity;

        for (var index = 0; index < FrequencyHz.Length; index++)
        {
            var frequency = FrequencyHz[index];
            if (frequency < left || frequency > right)
                continue;

            var power = PowerDb[index];
            if (power <= best_power)
                continue;

            best_power = power;
            best_index = index;
        }

        if (best_index < 0)
            return (0, double.NegativeInfinity);

        return (FrequencyHz[best_index], PowerDb[best_index]);
    }

    private static ComplexN[] FastFft(ComplexN[] Data)
    {
        var n = Data.Length;
        var result = new ComplexN[n];
        Array.Copy(Data, result, n);

        var j = 0;
        for (var i = 1; i < n; i++)
        {
            var bit = n >> 1;
            while ((j & bit) != 0)
            {
                j &= ~bit;
                bit >>= 1;
            }

            j |= bit;

            if (i < j)
                (result[i], result[j]) = (result[j], result[i]);
        }

        for (var len = 2; len <= n; len <<= 1)
        {
            var angle = -2 * Math.PI / len;
            var w_len = new ComplexN(Math.Cos(angle), Math.Sin(angle));

            for (var i = 0; i < n; i += len)
            {
                var w = ComplexN.One;
                var half_len = len / 2;
                for (var k = 0; k < half_len; k++)
                {
                    var u = result[i + k];
                    var v = result[i + k + half_len] * w;

                    result[i + k] = u + v;
                    result[i + k + half_len] = u - v;

                    w *= w_len;
                }
            }
        }

        return result;
    }
}
