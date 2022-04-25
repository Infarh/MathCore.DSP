namespace MathCore.DSP.Infrastructure;

internal class Tests
{
    public static double[] PolynomialCoefficients(int N)
    {
        var m = 1d;
        for (var i = 1; i <= N; i++) 
            m *= (N + i) / 2d;

        var a = new double[N + 1];
        (a[0], a[N]) = (m, 1);

        for (var i = 1; i < N; i++) 
            a[i] = a[i - 1] * 2 * (N - i + 1) / (2 * N - i + 1) / i;

        return a;
    }
}
