namespace MathCore.DSP.WindowFunctions;

public static class GaussianWindow
{
    public static double Value(int n, int N, double alpha) => Value((double)n / N, alpha);

    public static double Value(double x, double alpha) => Math.Exp(-2 * ((x - 0.5) / alpha).Pow2());
}
