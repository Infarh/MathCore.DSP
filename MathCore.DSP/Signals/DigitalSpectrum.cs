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
        if(index == index_double) return _Samples[index];
        var s1 = _Samples[index];
        var s2 = _Samples[index + 1];
        var k = index_double - index;
        return (1 - k) * s1 + k * s2;
    }
}