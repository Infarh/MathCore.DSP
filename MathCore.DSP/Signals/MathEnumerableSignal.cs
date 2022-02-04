namespace MathCore.DSP.Signals;

public static class MathEnumerableSignal
{
    private static IEnumerable<double> GetSamples(Func<double, double> f, double dt, int SamplesCount)
    {
        var x = 0d;
        var infinity = SamplesCount < 0;
        while (infinity || SamplesCount-- > 0)
            yield return f(x += dt);
    }

    public static EnumerableSignal Sin(double dt, double f0, int SamplesCount, double A = 1, double phi0 = 0)
    {
        var w0 = 2 * Math.PI * f0;
        return new EnumerableSignal(dt, GetSamples(t => A * Math.Sin(w0 * t + phi0), dt, SamplesCount));
    }

    public static EnumerableSignal Cos(double dt, double f0, int SamplesCount, double A = 1, double phi0 = 0)
    {
        var w0 = 2 * Math.PI * f0;
        return new EnumerableSignal(dt, GetSamples(t => A * Math.Sin(w0 * t + phi0), dt, SamplesCount));
    }

    public static EnumerableSignal Const(double dt, int SamplesCount, double A = 1) => new(dt, GetSamples(_ => A, dt, SamplesCount));
}