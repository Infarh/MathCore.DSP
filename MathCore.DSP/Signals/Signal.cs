namespace MathCore.DSP.Signals;

public static partial class Signal
{
    public static DigitalSignal Rect(double A, double t0, double t1, double t2, double dt, int SamplesCount, double A0 = 0)
    {
        var samples = new double[SamplesCount];

        for (var i = 0; i < SamplesCount; i++)
        {
            var t = t0 + i * dt;
            samples[i] = t >= t1 && t <= t2 ? A : A0;
        }

        return new SamplesDigitalSignal(dt, samples);
    }

    private static IEnumerable<double> Repeat(double Value, int Count)
    {
        for (var i = 0; i < Count; i++)
            yield return Value;
    }

    public static DigitalSignal Zero(double dt, double t0, int SamplesCount) => Const(0, dt, t0, SamplesCount);

    public static DigitalSignal Const(double A, double dt, double t0, int SamplesCount) => new EnumerableSignal(dt, Repeat(A, SamplesCount), t0);
}
