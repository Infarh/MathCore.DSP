namespace MathCore.DSP.WindowFunctions;
public class HammingWindow
{
    public static double Value(int n, int N, double a0) => Value((double)n / N, a0);

    public static double Value(double x, double a0) => a0 - (1 - a0) * Math.Cos(Consts.pi2 * x);

    private const double a0 = 25d / 46;
    private const double a1 = 1 - a0;

    public static double Value(int n, int N) => Value((double)n / N);

    public static double Value(double x) => a0 - a1 * Math.Cos(Consts.pi2 * x);
}
