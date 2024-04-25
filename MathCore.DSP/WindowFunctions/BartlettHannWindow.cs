namespace MathCore.DSP.WindowFunctions;

public static class BartlettHannWindow
{
    public static double Value(int n, int N) => Value((double)n / N);

    public static double Value(double x) => 0.62 - 0.48 * Math.Abs(x - 0.5) - 0.38 * Math.Cos(Consts.pi2 * x);
}
