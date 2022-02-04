using BenchmarkDotNet.Attributes;

using MathCore.DSP;

namespace Benchmarks;

[MemoryDiagnoser]
public class TestSum
{
    private static readonly int[] __Values = new int[1000000000];

    private static int SumFor(int[] values)
    {
        var sum = 0;
        for (var i = 0; i < values.Length; i++)
            sum += values[i];
        return sum;
    }

    private static int SumForOptimized(int[] values)
    {
        var sum = 0;
        for (var (i, count) = (0, values.Length); i < count; i++)
            sum += values[i];
        return sum;
    }

    private static int SumForeach(int[] values)
    {
        var sum = 0;
        foreach (var value in values)
            sum += value;
        return sum;
    }

    private static int SumLinq(int[] values)
    {
        return values.Sum();
    }

    private static int SumGoto(int[] values)
    {
        var sum = 0;
        var count = values.Length;
        var i = 0;
        goto Check;
    Iteration:
        sum += values[i];
        i++;
    Check:
        if (i < count)
            goto Iteration;

        return sum;
    }

    private static readonly Func<int[], int> __IlSummator = Summator.CreateDelegateInt();

    [Benchmark(Baseline = true)]
    public int SumFor()
    {
        return SumFor(__Values);
    }

    [Benchmark]
    public int SumForOptimized()
    {
        return SumForOptimized(__Values);
    }

    [Benchmark]
    public int SumForeach()
    {
        return SumForeach(__Values);
    }

    [Benchmark]
    public int SumLinq()
    {
        return SumLinq(__Values);
    }

    [Benchmark]
    public int SumGoto()
    {
        return SumGoto(__Values);
    }

    [Benchmark]
    public int SumIL()
    {
        return __IlSummator(__Values);
    }
}