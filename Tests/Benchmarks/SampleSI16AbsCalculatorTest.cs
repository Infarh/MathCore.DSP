using BenchmarkDotNet.Attributes;

using MathCore.DSP.Samples;

namespace Benchmarks;

public class SampleSI16AbsCalculatorTest
{
    private static readonly SampleSI16MagnitudeCalculator __Calculator = new();

    [Benchmark(Baseline = true)]
    public double SimpleCalculation()
    {
        var result = 0d;
        for (var i = -128; i <= 127; i++)
            for (var q = -128; q <= 127; q++)
            {
                var abs = MathF.Sqrt(i * (sbyte)i + q * (sbyte)q);
                result += abs;
            }

        return result;
    }

    [Benchmark]
    public double Calculator()
    {
        var result = 0d;
        for (var i = -128; i <= 127; i++)
            for (var q = -128; q <= 127; q++)
            {
                var abs = __Calculator.GetMagnitude(new((sbyte)i, (sbyte)q));
                result += abs;
            }

        return result;
    }
}