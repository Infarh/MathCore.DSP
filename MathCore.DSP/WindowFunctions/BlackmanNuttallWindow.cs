namespace MathCore.DSP.WindowFunctions;

public static class BlackmanNuttallWindow
{
    private const double a0 = 0.3635819;
    private const double a1 = 0.4891775;
    private const double a2 = 0.1365995;
    private const double a3 = 0.0106411;

    private const double pi4 = Consts.pi2 * 2;
    private const double pi6 = Consts.pi2 * 3;

    public static double Value(int n, int N) => Value((double)n / N);

    public static double Value(double x) => 
        a0
        - a1 * Math.Cos(Consts.pi2 * x)
        + a2 * Math.Cos(Consts.pi2 * 2 * x)
        - a3 * Math.Cos(Consts.pi2 * 3 * x);

    public static double Value(int n, int N, double a0, double a1, double a2, double a3) => Value((double)n / N, a0, a1, a2, a3);

    public static double Value(double x, double a0, double a1, double a2, double a3) => 
        a0
        - a1 * Math.Cos(Consts.pi2 * x)
        + a2 * Math.Cos(pi4 * x)
        - a3 * Math.Cos(pi6 * x);
}
