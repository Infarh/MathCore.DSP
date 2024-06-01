namespace MathCore.DSP.WindowFunctions;

public static class KaiserWindow
{
    public static double Value(int n, int N, double alpha) => Value((double)n / N, alpha);

    public static double Value(double x, double alpha) => 
        SpecialFunctions.Bessel.I0(alpha * (1 - (2 * (x - 0.5)).Pow2()).Sqrt()) /
        SpecialFunctions.Bessel.I0(alpha);
}
