namespace MathCore.DSP.Filters.Builders;

public readonly record struct HighPassBuilder(double dt)
{
    public HighPassButterworthBuilder Butterworth(double fs, double fp, double Gp, double Gs) => new(dt, fs, fp, Gp, Gs);

    public HighPassChebyshevBuilder Chebyshev(double fs, double fp, double Gp, double Gs) => new(dt, fs, fp, Gp, Gs);
}

public readonly record struct HighPassButterworthBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    public ButterworthHighPass Create() => new(dt, fp, fs, Gp, Gs);

    public static implicit operator Filter(HighPassButterworthBuilder builder) => builder.Create();
}

public readonly record struct HighPassChebyshevBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    public ChebyshevHighPass Create() => new(dt, fp, fs, Gp, Gs);

    public static implicit operator Filter(HighPassChebyshevBuilder builder) => builder.Create();
}

public readonly record struct HighPassEllipticBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    public EllipticHighPass Create() => new(dt, fp, fs, Gp, Gs);

    public static implicit operator Filter(HighPassEllipticBuilder builder) => builder.Create();
}
