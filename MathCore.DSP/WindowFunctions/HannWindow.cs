namespace MathCore.DSP.WindowFunctions;

public static class HannWindow
{
    public static double Value(int n, int N) => Value((double)n / N);

    public static double Value(double x) => Math.Sin(Math.PI * x).Pow2();
    //public static double Value(int n, int N) => 0.5 - 0.5 * Math.Cos(Math.PI * n / N);
}
