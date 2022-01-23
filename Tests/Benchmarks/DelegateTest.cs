using BenchmarkDotNet.Attributes;

namespace Benchmarks;

[MemoryDiagnoser]
public class DelegateTest
{
    private static int _X = 5;
    private static int _Y = 7;

    private static int Sum(int a, int b) => a + b;

    private static readonly Func<int, int, int> __FuncSum = Sum;

    private static readonly Delegate __DelegateSum = new Func<int, int, int>(Sum);

    [Benchmark(Baseline = true)]
    public int FuncSum()
    {
        var result = 0;
        var sum = __FuncSum;
        for (var i = 0; i < 10000000; i++)
            result = sum(_X, _Y);
        return result;
    }

    [Benchmark]
    public int DelegateOptimizedSum()
    {
        var result = 0;
        var sum = (Func<int, int, int>)__DelegateSum;
        for (var i = 0; i < 10000000; i++)
            result = sum(_X, _Y)!;
        return result;
    }
}