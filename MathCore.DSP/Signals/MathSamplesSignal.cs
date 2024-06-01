namespace MathCore.DSP.Signals;

public static class MathSamplesSignal
{
    public static SamplesDigitalSignal Sin(double A, double f0, double phi0, double dt, int SamplesCount)
    {
        var w0      = Consts.pi2 * f0 * dt;
        var samples = new double[SamplesCount];
        for (var i = 0; i < SamplesCount; i++)
            samples[i] = A * Math.Sin(w0 * i + phi0);
        return new(dt, samples);
    }

    public static SamplesDigitalSignal Cos(double A, double f0, double phi0, double dt, int SamplesCount)
    {
        var w0      = Consts.pi2 * f0 * dt;
        var samples = new double[SamplesCount];
        for (var i = 0; i < SamplesCount; i++)
            samples[i] = A * Math.Cos(w0 * i + phi0);
        return new(dt, samples);
    }
}