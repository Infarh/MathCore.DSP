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

    [TestMethod]
    public void FilterIIRInterleaved_MatchesSampleSI16Filtering()
    {
        var random = new Random(17);
        var raw_iq = new byte[256 * 2];
        random.NextBytes(raw_iq);

        var samples = new SampleSI16[raw_iq.Length / 2];
        for (var i = 0; i < samples.Length; i++)
            samples[i] = new(unchecked((sbyte)raw_iq[2 * i]), unchecked((sbyte)raw_iq[2 * i + 1]));

        var a = new[] { 1.0, -0.78, 0.16 };
        var b = new[] { 0.15, 0.05, 0.01 };

        var expected_samples = samples.AsEnumerable().FilterIIR(a, b).ToArray();

        var actual_raw = raw_iq.AsSpan().FilterIIRInterleaved(a, b);

        Assert.AreEqual(raw_iq.Length, actual_raw.Length);
        for (var i = 0; i < expected_samples.Length; i++)
        {
            Assert.AreEqual(unchecked((byte)expected_samples[i].I), actual_raw[2 * i], $"Неверная I-компонента на позиции {i}");
            Assert.AreEqual(unchecked((byte)expected_samples[i].Q), actual_raw[2 * i + 1], $"Неверная Q-компонента на позиции {i}");
        }
    }

    [TestMethod]
    public void FilterIIRInterleaved_InPlaceAndOutPlace_ReturnSameResult()
    {
        var random = new Random(23);
        var raw_iq = new byte[128 * 2];
        random.NextBytes(raw_iq);

        var a = new[] { 1.0, -0.62 };
        var b = new[] { 0.25, 0.11 };

        var expected = raw_iq.AsSpan().FilterIIRInterleaved(a, b);

        var actual = raw_iq.ToArray();
        actual.AsSpan().FilterIIRInterleaved(a, b, new double[a.Length], new double[a.Length]);

        CollectionAssert.AreEqual(expected, actual);
    }

    [TestMethod]
    public void FilterIIRInterleaved_ThrowsIfInputLengthIsOdd()
    {
        var raw_iq = new byte[] { 1, 2, 3 };
        var destination = new byte[raw_iq.Length];
        var a = new[] { 1.0, -0.5 };
        var b = new[] { 0.5, 0.1 };

        try
        {
            raw_iq.AsSpan().FilterIIRInterleaved(destination, a, b, new double[a.Length], new double[a.Length]);
            Assert.Fail("Ожидалось исключение InvalidOperationException");
        }
        catch (InvalidOperationException)
        {
            // Ожидаемое поведение
        }
    }

    [TestMethod]
    public void FilterIIRInterleaved_Stream_MatchesSpanProcessing()
    {
        var random = new Random(31);
        var raw_iq = new byte[257 * 2];
        random.NextBytes(raw_iq);

        var a = new[] { 1.0, -0.81, 0.21 };
        var b = new[] { 0.11, 0.05, 0.02 };

        var expected = raw_iq.AsSpan().FilterIIRInterleaved(a, b);

        using var source = new MemoryStream(raw_iq, writable: false);
        using var destination = new MemoryStream();

        var written = source.FilterIIRInterleaved(destination, new byte[15], a, b, new double[a.Length], new double[a.Length]);
        var actual = destination.ToArray();

        Assert.AreEqual(expected.Length, (int)written);
        CollectionAssert.AreEqual(expected, actual);
    }

    [TestMethod]
    public void FilterIIRInterleaved_Stream_ThrowsIfTotalLengthIsOdd()
    {
        var raw_iq = new byte[] { 1, 2, 3, 4, 5 };
        var a = new[] { 1.0, -0.5 };
        var b = new[] { 0.5, 0.1 };

        using var source = new MemoryStream(raw_iq, writable: false);
        using var destination = new MemoryStream();

        try
        {
            source.FilterIIRInterleaved(destination, new byte[4], a, b, new double[a.Length], new double[a.Length]);
            Assert.Fail("Ожидалось исключение IOException");
        }
        catch (IOException)
        {
            // Ожидаемое поведение
        }
    }

    [TestMethod]
    public async Task FilterIIRInterleavedAsync_Stream_MatchesSyncVersion()
    {
        var random = new Random(37);
        var raw_iq = new byte[333 * 2];
        random.NextBytes(raw_iq);

        var a = new[] { 1.0, -0.74, 0.19 };
        var b = new[] { 0.14, 0.06, 0.01 };

        using var source_sync = new MemoryStream(raw_iq, writable: false);
        using var destination_sync = new MemoryStream();
        source_sync.FilterIIRInterleaved(destination_sync, new byte[13], a, b, new double[a.Length], new double[a.Length]);
        var expected = destination_sync.ToArray();

        using var source_async = new MemoryStream(raw_iq, writable: false);
        using var destination_async = new MemoryStream();
        var written = await source_async.FilterIIRInterleavedAsync(destination_async, new byte[13], a, b, new double[a.Length], new double[a.Length]);
        var actual = destination_async.ToArray();

        Assert.AreEqual(expected.Length, (int)written);
        CollectionAssert.AreEqual(expected, actual);
    }

    [TestMethod]
    public async Task FilterIIRInterleavedAsync_Stream_ThrowsIfTotalLengthIsOdd()
    {
        var raw_iq = new byte[] { 1, 2, 3, 4, 5 };
        var a = new[] { 1.0, -0.5 };
        var b = new[] { 0.5, 0.1 };

        using var source = new MemoryStream(raw_iq, writable: false);
        using var destination = new MemoryStream();

        try
        {
            await source.FilterIIRInterleavedAsync(destination, new byte[4], a, b, new double[a.Length], new double[a.Length]);
            Assert.Fail("Ожидалось исключение IOException");
        }
        catch (IOException)
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