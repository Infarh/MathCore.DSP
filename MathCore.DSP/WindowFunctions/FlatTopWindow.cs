namespace MathCore.DSP.WindowFunctions;

public static class FlatTopWindow
{
    private const double a1 = 1.93;
    private const double a2 = 1.29;
    private const double a3 = 0.388;
    private const double a4 = 0.032;

    private const double pi2 = Consts.pi2 * 1;
    private const double pi4 = Consts.pi2 * 2;
    private const double pi6 = Consts.pi2 * 3;
    private const double pi8 = Consts.pi2 * 4;

    public static double Value(int n, int N) => Value((double)n / N);

    public static double Value(double x) =>
        1
        - a1 * Math.Cos(pi2 * x)
        + a2 * Math.Cos(pi4 * x)
        - a3 * Math.Cos(pi6 * x)
        + a4 * Math.Cos(pi8 * x);

    public static double Value(int n, int N, double a1, double a2, double a3, double a4) => Value((double)n / N, a1, a2, a3, a4);

    public static double Value(double x, double a1, double a2, double a3, double a4) =>
        1
        - a1 * Math.Cos(pi2 * x)
        + a2 * Math.Cos(pi4 * x)
        - a3 * Math.Cos(pi6 * x)
        + a4 * Math.Cos(pi8 * x);
}
