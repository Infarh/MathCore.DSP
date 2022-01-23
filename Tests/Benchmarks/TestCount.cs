using System.Collections;

using BenchmarkDotNet.Attributes;
// ReSharper disable ArrayCount
// ReSharper disable IListCount

namespace Benchmarks;

[MemoryDiagnoser]
public class TestCount
{
    private static readonly int[] __Values = new int[1000000000];

    private static int GetLength(int[] values) => values.Length;

    private static int GetCount(int[] values) => values.Count();

    private static int GetListTCount(IList<int> values) => values.Count;
    private static int GetListCount(IList values) => values.Count;
    private static int GetCollectionTCount(ICollection<int> values) => values.Count;
    private static int GetCollectionCount(ICollection<int> values) => values.Count;

    [Benchmark(Baseline = true)]
    public int Length()
    {
        var count = 0;
        for (var i = 0; i < 1000000; i++)
            count = GetLength(__Values);
        return count;
    }

    [Benchmark]
    public int Count()
    {
        var count = 0;
        for (var i = 0; i < 1000000; i++)
            count = GetCount(__Values);
        return count;
    }

    [Benchmark]
    public int ListTCount()
    {
        var count = 0;
        for (var i = 0; i < 1000000; i++)
            count = GetListTCount(__Values);
        return count;
    }

    [Benchmark]
    public int ListCount()
    {
        var count = 0;
        for (var i = 0; i < 1000000; i++)
            count = GetListCount(__Values);
        return count;
    }

    [Benchmark]
    public int CollectionTCount()
    {
        var count = 0;
        for (var i = 0; i < 1000000; i++)
            count = GetCollectionTCount(__Values);
        return count;
    }

    [Benchmark]
    public int CollectionCount()
    {
        var count = 0;
        for (var i = 0; i < 1000000; i++)
            count = GetCollectionCount(__Values);
        return count;
    }
}