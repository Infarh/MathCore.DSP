using System.Data;

using MathCore.DSP.Fourier;

namespace MathCore.DSP.Signals;

public class DigitalSpectrum
{
    private double _fd;
    private readonly Complex[] _Samples;
    private readonly double _t0;

    public double fd => _fd;

    public double df => _fd / _Samples.Length;

    public double t0 => _t0;

    public int SamplesCount => _Samples.Length;

    public ref Complex this[int m] => ref _Samples[m];

    public DigitalSpectrum(double fd, Complex[] Samples, double t0 = 0)
    {
        _fd = fd;
        _Samples = Samples;
        _t0 = t0;
    }

    public static DigitalSpectrum operator /(DigitalSpectrum X, DigitalSpectrum Y)
    {
        var x_samples = X._Samples;
        var y_samples = Y._Samples;

        var result_samples = new Complex[x_samples.Length];

        for (var i = 0; i < result_samples.Length; i++)
            result_samples[i] = x_samples[i] / y_samples[i];

        return new(X._fd, result_samples, X._t0);
    }

    public Complex GetValue(double f)
    {
        if (f > fd || f < 0) throw new ArgumentOutOfRangeException(nameof(f), f, "Частота должна быть в диапазоне от 0 до fd");

        var index_double = f / df;
        var index = (int)index_double;
        if (index == index_double) return _Samples[index];
        var s1 = _Samples[index];
        var s2 = _Samples[index + 1];
        var k = index_double - index;
        return (1 - k) * s1 + k * s2;
    }

    public DigitalSignal GetSignalReal(bool ThrowIfImExists = false, double ImPowerThreshold_dB = 40)
    {
        var complex_samples = _Samples.FastFourierInverse();
        var (re, im) = complex_samples.ToReIm();

        if (ThrowIfImExists)
        {
            var re_power = re.Average(v => v * v);
            var im_power = im.Average(v => v * v);

            var im_power_db = im_power.In_dB_byPower();
            var re_power_db = re_power.In_dB_byPower();
            var im_to_re_power_db = re_power_db - im_power_db;

            if (im_to_re_power_db < ImPowerThreshold_dB)
                throw new InvalidOperationException($"В процессе выполнения обратного преобразования Фурье от спектра была получена комплексная составляющая, мощность которой превысила порог в {ImPowerThreshold_dB} дБ")
                {
                    Data =
                    {
                        { "ImPower", im_power},
                        { "ImPower dB", im_power_db},
                        { "RePower", re_power},
                        { "RePower dB", re_power_db},
                        { "RePower dB", re_power_db},
                        { "Im / Re Power dB", im_to_re_power_db},
                        { "ImPowerThreshold dB", ImPowerThreshold_dB},
                    }
                };
        }

        return new SamplesDigitalSignal(1 / fd, re, _t0);
    }

    public (DigitalSignal I, DigitalSignal Q) GetSignalComplex()
    {
        var complex_samples = _Samples.FastFourierInverse();
        var (re, im) = complex_samples.ToReIm();

        var dt = 1 / _fd;
        var i = new SamplesDigitalSignal(dt, re, _t0);
        var q = new SamplesDigitalSignal(dt, im, _t0);
        return (i, q);
    }
}