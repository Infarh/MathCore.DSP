namespace MathCore.DSP.Filters.Builders;

public readonly record struct LowPassBuilder(double dt)
{
    public LowPassButterworthBuilder Butterworth(double fs, double fp, double Gp, double Gs) => new(dt, fs, fp, Gp, Gs);

    public LowPassChebyshevBuilder Chebyshev(double fs, double fp, double Gp, double Gs) => new(dt, fs, fp, Gp, Gs);

    public LowPassRCBuilder RC(double f0) => new(dt, f0);

    public LowPassRCExponentialBuilder RCExponential(double f0) => new(dt, f0);
}

public readonly record struct LowPassButterworthBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    public ButterworthLowPass Create() => new(dt, fp, fs, Gp, Gs);

    public static implicit operator Filter(LowPassButterworthBuilder builder) => builder.Create();
}

public readonly record struct LowPassChebyshevBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    public ChebyshevLowPass Create() => new(dt, fp, fs, Gp, Gs);

    public static implicit operator Filter(LowPassChebyshevBuilder builder) => builder.Create();
}

public readonly record struct LowPassEllipticBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    public EllipticLowPass Create() => new(dt, fp, fs, Gp, Gs);

    public static implicit operator Filter(LowPassEllipticBuilder builder) => builder.Create();
}

public readonly record struct LowPassRCBuilder(double dt, double f0)
{
    public RCLowPass Create() => new(dt, f0);

    public static implicit operator Filter(LowPassRCBuilder builder) => builder.Create();
}

public readonly record struct LowPassRCExponentialBuilder(double dt, double f0)
{
    public RCExponentialLowPass Create() => new(dt, f0);

    public static implicit operator Filter(LowPassRCExponentialBuilder builder) => builder.Create();
}

