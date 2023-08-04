namespace MathCore.DSP.Tests.Infrastructure;

public static class AccuracyComplex
{
    public static IEqualityComparer<Complex> Eps(double eps) => new AccuracyComplexEqualityComparer(eps);
}

public readonly struct AccuracyComplexEqualityComparer :
    IEqualityComparer<Complex>
{
    private double Eps { get; init; }

    public AccuracyComplexEqualityComparer(double Eps)
    {
        if (Eps < 0.0)
            throw new ArgumentOutOfRangeException(nameof(Eps), (object)Eps, "Значение точности не должно быть меньше нуля");
        this.Eps = !double.IsNaN(Eps) ? Eps : throw new ArgumentException("Значение точности не должно быть NaN", nameof(Eps));
    }

    public bool Equals(Complex x, Complex y)
    {
        var eps = Eps;
        var (delta_re, delta_im) = x - y;
        return Math.Abs(delta_re) <= eps 
            && Math.Abs(delta_im) <= eps;
    }

    public int GetHashCode(Complex z)
    {
        if (z.IsNaN) return z.GetHashCode();
        var (re, im) = z;
        var eps = Eps;
        var value = new Complex(
            Math.Round(re * eps) / eps,
            Math.Round(im * eps) / eps
        );
        return value.GetHashCode();
    }
}


