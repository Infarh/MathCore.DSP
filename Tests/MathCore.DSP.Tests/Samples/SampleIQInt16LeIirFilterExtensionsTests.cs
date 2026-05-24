using System.Buffers.Binary;
using System.IO;

using MathCore.DSP.Samples.Extensions;

namespace MathCore.DSP.Tests.Samples;

[TestClass]
public class SampleIQInt16LeIirFilterExtensionsTests
{
    [TestMethod]
    public void FilterIIRInterleavedInt16Le_MatchesIndependentChannelFiltering()
    {
        var random = new Random(107);

        var samples_count = 256;
        var raw = new byte[samples_count * 4];
        var i_values = new short[samples_count];
        var q_values = new short[samples_count];

        for (var i = 0; i < samples_count; i++)
        {
            i_values[i] = (short)random.Next(short.MinValue, short.MaxValue + 1);
            q_values[i] = (short)random.Next(short.MinValue, short.MaxValue + 1);

            BinaryPrimitives.WriteInt16LittleEndian(raw.AsSpan(i * 4, 2), i_values[i]);
            BinaryPrimitives.WriteInt16LittleEndian(raw.AsSpan(i * 4 + 2, 2), q_values[i]);
        }

        var a = new[] { 1.0, -0.85, 0.2 };
        var b = new[] { 0.12, 0.06, 0.01 };

        var expected_i = i_values.Select(v => (double)v).FilterIIR(a, b).ToArray();
        var expected_q = q_values.Select(v => (double)v).FilterIIR(a, b).ToArray();

        var actual_raw = raw.AsSpan().FilterIIRInterleavedInt16Le(a, b);

        for (var i = 0; i < samples_count; i++)
        {
            var actual_i = BinaryPrimitives.ReadInt16LittleEndian(actual_raw.AsSpan(i * 4, 2));
            var actual_q = BinaryPrimitives.ReadInt16LittleEndian(actual_raw.AsSpan(i * 4 + 2, 2));

            Assert.AreEqual(ClampToInt16(expected_i[i]), actual_i, $"Неверная I-компонента на позиции {i}");
            Assert.AreEqual(ClampToInt16(expected_q[i]), actual_q, $"Неверная Q-компонента на позиции {i}");
        }
    }

    [TestMethod]
    public void FilterIIRInterleavedInt16Le_ThrowsIfLengthNotMultipleOf4()
    {
        var raw = new byte[] { 1, 2, 3, 4, 5, 6 };
        var destination = new byte[raw.Length];

        var a = new[] { 1.0, -0.6 };
        var b = new[] { 0.2, 0.1 };

        try
        {
            raw.AsSpan().FilterIIRInterleavedInt16Le(destination, a, b, new double[a.Length], new double[a.Length]);
            Assert.Fail("Ожидалось исключение InvalidOperationException");
        }
        catch (InvalidOperationException)
        {
            // Ожидаемое поведение
        }
    }

    [TestMethod]
    public void FilterIIRInterleavedInt16Le_Stream_MatchesSpanProcessing()
    {
        var random = new Random(109);
        var raw = new byte[255 * 4];
        random.NextBytes(raw);

        var a = new[] { 1.0, -0.81, 0.18 };
        var b = new[] { 0.13, 0.05, 0.01 };

        var expected = raw.AsSpan().FilterIIRInterleavedInt16Le(a, b);

        using var source = new MemoryStream(raw, writable: false);
        using var destination = new MemoryStream();

        var written = source.FilterIIRInterleavedInt16Le(destination, new byte[11], a, b, new double[a.Length], new double[a.Length]);
        var actual = destination.ToArray();

        Assert.AreEqual(expected.Length, (int)written);
        CollectionAssert.AreEqual(expected, actual);
    }

    [TestMethod]
    public async Task FilterIIRInterleavedInt16LeAsync_Stream_MatchesSyncVersion()
    {
        var random = new Random(113);
        var raw = new byte[199 * 4];
        random.NextBytes(raw);

        var a = new[] { 1.0, -0.77, 0.16 };
        var b = new[] { 0.12, 0.06, 0.02 };

        using var source_sync = new MemoryStream(raw, writable: false);
        using var destination_sync = new MemoryStream();
        source_sync.FilterIIRInterleavedInt16Le(destination_sync, new byte[9], a, b, new double[a.Length], new double[a.Length]);
        var expected = destination_sync.ToArray();

        using var source_async = new MemoryStream(raw, writable: false);
        using var destination_async = new MemoryStream();
        var written = await source_async.FilterIIRInterleavedInt16LeAsync(destination_async, new byte[9], a, b, new double[a.Length], new double[a.Length]);
        var actual = destination_async.ToArray();

        Assert.AreEqual(expected.Length, (int)written);
        CollectionAssert.AreEqual(expected, actual);
    }

    private static short ClampToInt16(double value) => value switch
    {
        > short.MaxValue => short.MaxValue,
        < short.MinValue => short.MinValue,
        _ => (short)Math.Round(value)
    };
}
