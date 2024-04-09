using static System.Math;


namespace MathCore.DSP.WindowFunctions;

// https://habr.com/ru/post/514170/
public class BlackmanWindow
{
    private const double alpha = 0.16;

    private const double a0 = (1 - alpha) / 2;
    private const double a1 = alpha / 2;

    private const double pi4 = Consts.pi2 * 2;

    public static double Value(int n, int N, double a) => Value((double)n / N, a);

    public static double Value(double x, double a) => (1 - a) - Cos(Consts.pi2 * x) + a * Cos(pi4 * x);

    public static double Value(int n, int N) => Value((double)n / N);

    public static double Value(double x) => a0 - 0.5 * Cos(Consts.pi2 * x) + a1 * Cos(pi4 * x);

    //public static double Value(double x) =>
    //    Abs(x) <= 0.5
    //        ? (25 * Cos(2 * PI * x) + 4 * Cos(4 * PI * x) + 21) / 50
    //        : 0;

    //public static double Nuttall(double x) =>
    //    Abs(x) <= 0.5
    //        ? (121_849 * Cos(2 * PI * x) + 36_058 * Cos(4 * PI * x) + 3_151 * Cos(6 * PI * x) + 88_942) / 250_000
    //        : 0;
}