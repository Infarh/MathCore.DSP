using MathCore.DSP.Samples;
using MathCore.DSP.Samples.Extensions;

namespace MathCore.DSP.Tests.Samples;

[TestClass]
public class SampleSI16IirFilterExtensionsTests
{
    [TestMethod]
    public void FilterIIR_Enumerable_MatchesIndependentDoubleFilteringForIAndQ()
    {
        var random = new Random(7);
        var samples = Enumerable
           .Range(0, 256)
           .Select(_ => new SampleSI16((sbyte)random.Next(sbyte.MinValue, sbyte.MaxValue + 1), (sbyte)random.Next(sbyte.MinValue, sbyte.MaxValue + 1)))
           .ToArray();

        var a = new[] { 1.0, -0.82, 0.18 };
        var b = new[] { 0.12, 0.07, 0.02 };

        var expected_i = samples.Select(s => (double)s.I).FilterIIR(a, b).ToArray();
        var expected_q = samples.Select(s => (double)s.Q).FilterIIR(a, b).ToArray();

        var actual = samples.AsEnumerable().FilterIIR(a, b).ToArray();

        Assert.AreEqual(samples.Length, actual.Length);

        for (var i = 0; i < actual.Length; i++)
        {
            var expected_i_sample = ClampToSByte(expected_i[i]);
            var expected_q_sample = ClampToSByte(expected_q[i]);

            Assert.AreEqual(expected_i_sample, actual[i].I, $"Неверная I-компонента на позиции {i}");
            Assert.AreEqual(expected_q_sample, actual[i].Q, $"Неверная Q-компонента на позиции {i}");
        }
    }

    [TestMethod]
    public void FilterIIR_SpanAndEnumerable_ReturnSameResult()
    {
        var random = new Random(11);
        var samples = Enumerable
           .Range(0, 128)
           .Select(_ => new SampleSI16((sbyte)random.Next(sbyte.MinValue, sbyte.MaxValue + 1), (sbyte)random.Next(sbyte.MinValue, sbyte.MaxValue + 1)))
           .ToArray();

        var a = new[] { 1.0, -0.7 };
        var b = new[] { 0.3, 0.15 };

        var expected = samples.AsEnumerable().FilterIIR(a, b).ToArray();

        var actual = new SampleSI16[samples.Length];
        samples.AsSpan().FilterIIR(actual, a, b, new double[a.Length], new double[a.Length]);

        CollectionAssert.AreEqual(expected, actual);
    }

    [TestMethod]
    public void FilterIIR_ThrowsIfDestinationShorterThanInput()
    {
        var samples = new SampleSI16[]
        {
            new(1, 2),
            new(3, 4),
            new(5, 6)
        };

        var destination = new SampleSI16[2];
        var a = new[] { 1.0, -0.5 };
        var b = new[] { 0.5, 0.1 };

        try
        {
            samples.AsSpan().FilterIIR(destination, a, b, new double[a.Length], new double[a.Length]);
            Assert.Fail("Ожидалось исключение ArgumentException");
        }
        catch (ArgumentException)
        {
            // Ожидаемое поведение
        }
    }

    private static sbyte ClampToSByte(double value) => value switch
    {
        > sbyte.MaxValue => sbyte.MaxValue,
        < sbyte.MinValue => sbyte.MinValue,
        _ => (sbyte)Math.Round(value)
    };
}