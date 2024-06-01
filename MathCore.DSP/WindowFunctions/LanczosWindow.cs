namespace MathCore.DSP.WindowFunctions;

/// <summary>Окно Ланцоша</summary>
public static class LanczosWindow
{
    public static double Value(int n, int N) => Value((double)n / N);
    public static double Value(double x) => MathEx.Sinc(Consts.pi2 * x - Consts.pi);

    //public static double Value(int n, int N) => MathEx.Sinc(Consts.pi2 * n / N - Consts.pi);
}
