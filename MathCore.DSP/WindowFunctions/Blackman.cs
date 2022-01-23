using static System.Math;


namespace MathCore.DSP.WindowFunctions;

// https://habr.com/ru/post/514170/
public static class Functions
{
    public static double Blackman(double x) =>
        Abs(x) <= 0.5
            ? (25 * Cos(2 * PI * x) + 4 * Cos(4 * PI * x) + 21) / 50
            : 0;

    public static double Nuttall(double x) =>
        Abs(x) <= 0.5
            ? (121_849 * Cos(2 * PI * x) + 36_058 * Cos(4 * PI * x) + 3_151 * Cos(6 * PI * x) + 88_942) / 250_000
            : 0;


}