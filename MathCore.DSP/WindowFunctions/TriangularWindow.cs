namespace MathCore.DSP.WindowFunctions;

public static class TriangularWindow
{
    public static double Value(int n, int N) => Value((double)n / N);

    public static double Value(double x) => 1 - 2 * Math.Abs(x - 0.5);

    //public static double Value(int n, int N) => 1 - 2 * Math.Abs((double)n / N - 0.5) ;
    //public static double Value(int n, int N) => (N - 2 * Math.Abs(n - N / 2d)) / N;
}
