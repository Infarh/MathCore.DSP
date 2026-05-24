using System.IO;

using MathCore.DSP.Samples.Extensions;

namespace MathCore.DSP.Tests.Samples;

[TestClass]
public class SampleIQInt16LeFirFilterExtensionsTests
{
    [TestMethod]
    public void FilterFIRInterleavedInt16Le_MatchesIndependentChannelFiltering()
    {
        var random = new Random(127);

        var samples_count = 256;
        var raw = new byte[samples_count * 4];
        var i_values = new short[samples_count];
        var q_values = new short[samples_count];

        for (var i = 0; i < samples_count; i++)
        {
            i_values[i] = (short)random.Next(short.MinValue, short.MaxValue + 1);
            q_values[i] = (short)random.Next(short.MinValue, short.MaxValue + 1);

            System.Buffers.Binary.BinaryPrimitives.WriteInt16LittleEndian(raw.AsSpan(i * 4, 2), i_values[i]);
            System.Buffers.Binary.BinaryPrimitives.WriteInt16LittleEndian(raw.AsSpan(i * 4 + 2, 2), q_values[i]);
        }

        double[] impulse_response = [0.2, 0.3, 0.2, 0.1];

        var expected_i = i_values.Select(v => (double)v).FilterFIR(impulse_response).ToArray();
        var expected_q = q_values.Select(v => (double)v).FilterFIR(impulse_response).ToArray();

        var actual_raw = raw.AsSpan().FilterFIRInterleavedInt16Le(impulse_response);

        for (var i = 0; i < samples_count; i++)
        {
            var actual_i = System.Buffers.Binary.BinaryPrimitives.ReadInt16LittleEndian(actual_raw.AsSpan(i * 4, 2));
            var actual_q = System.Buffers.Binary.BinaryPrimitives.ReadInt16LittleEndian(actual_raw.AsSpan(i * 4 + 2, 2));

            Assert.AreEqual(ClampToInt16(expected_i[i]), actual_i, $"Неверная I-компонента на позиции {i}");
            Assert.AreEqual(ClampToInt16(expected_q[i]), actual_q, $"Неверная Q-компонента на позиции {i}");
        }
    }

    [TestMethod]
    public void FilterFIRInterleavedInt16Le_Stream_MatchesSpanProcessing()
    {
        var random = new Random(131);
        var raw = new byte[111 * 4];
        random.NextBytes(raw);

        double[] impulse_response = [0.12, 0.25, 0.2, 0.15, 0.1];

        var expected = raw.AsSpan().FilterFIRInterleavedInt16Le(impulse_response);

        using var source = new MemoryStream(raw, writable: false);
        using var destination = new MemoryStream();

        var written = source.FilterFIRInterleavedInt16Le(destination, new byte[7], impulse_response, new double[impulse_response.Length], new double[impulse_response.Length]);
        var actual = destination.ToArray();

        Assert.AreEqual(expected.Length, (int)written);
        CollectionAssert.AreEqual(expected, actual);
    }

    [TestMethod]
    public async Task FilterFIRInterleavedInt16LeAsync_Stream_MatchesSyncVersion()
    {
        var random = new Random(137);
        var raw = new byte[95 * 4];
        random.NextBytes(raw);

        double[] impulse_response = [0.18, 0.27, 0.21, 0.13];

        using var source_sync = new MemoryStream(raw, writable: false);
        using var destination_sync = new MemoryStream();
        source_sync.FilterFIRInterleavedInt16Le(destination_sync, new byte[9], impulse_response, new double[impulse_response.Length], new double[impulse_response.Length]);
        var expected = destination_sync.ToArray();

        using var source_async = new MemoryStream(raw, writable: false);
        using var destination_async = new MemoryStream();
        var written = await source_async.FilterFIRInterleavedInt16LeAsync(destination_async, new byte[9], impulse_response, new double[impulse_response.Length], new double[impulse_response.Length]);
        var actual = destination_async.ToArray();

        Assert.AreEqual(expected.Length, (int)written);
        CollectionAssert.AreEqual(expected, actual);
    }

    [TestMethod]
    public void FilterFIRInterleavedInt16Le_ThrowsIfLengthNotMultipleOf4()
    {
        var raw = new byte[] { 1, 2, 3, 4, 5, 6 };
        var destination = new byte[raw.Length];

        double[] impulse_response = [0.2, 0.1, 0.05];

        try
        {
            raw.AsSpan().FilterFIRInterleavedInt16Le(destination, impulse_response, new double[impulse_response.Length], new double[impulse_response.Length]);
            Assert.Fail("Ожидалось исключение InvalidOperationException");
        }
        catch (InvalidOperationException)
        {
            // Ожидаемое поведение
        }
    }

    private static short ClampToInt16(double value) => value switch
    {
        > short.MaxValue => short.MaxValue,
        < short.MinValue => short.MinValue,
        _ => (short)Math.Round(value)
    };
}
