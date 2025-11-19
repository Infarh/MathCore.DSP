using BenchmarkDotNet.Attributes;

using MathCore.DSP.Samples;

namespace Benchmarks;

public class SampleSI16ArgCalculatorTest
{
    private static readonly SampleSI16ArgumentCalculator __Calculator = new();

    [Benchmark(Baseline = false)]
    public double SimpleCalculation()
    {
        var result = 0d;
        for (var i = -128; i <= 127; i++)
            for (var q = -128; q <= 127; q++)
            {
                var abs = MathF.Atan2((sbyte)q, (sbyte)i);
                result += abs;
            }

        return result;
    }

    [Benchmark(Baseline = true)]
    public double Calculator()
    {
        var result = 0d;
        for (var i = -128; i <= 127; i++)
            for (var q = -128; q <= 127; q++)
            {
                var abs = __Calculator.GetArgument(new((sbyte)i, (sbyte)q));
                result += abs;
            }

        return result;
    }
}