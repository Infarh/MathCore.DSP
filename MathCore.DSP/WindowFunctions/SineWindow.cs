namespace MathCore.DSP.WindowFunctions;

public static class SineWindow
{
    public static double Value(int n, int N) => Value((double)n / N);

    public static double Value(double x) => Math.Sin(Math.PI * x);
}
