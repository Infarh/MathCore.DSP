using System.Buffers.Binary;
using System.IO;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;

namespace MathCore.DSP.Samples.Extensions;

/// <summary>Методы КИХ-фильтрации для interleaved IQ-потока int16 little-endian (USRP B210)</summary>
public static class SampleIQInt16LeFirFilterExtensions
{
    /// <summary>Отфильтровать raw буфер interleaved IQ int16 little-endian</summary>
    /// <param name="Samples">Входной буфер формата I16LE,Q16LE,I16LE,Q16LE</param>
    /// <param name="Destination">Выходной буфер того же формата</param>
    /// <param name="ImpulseResponse">Массив коэффициентов КИХ-фильтра</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    public static void FilterFIRInterleavedInt16Le(
        this ReadOnlySpan<byte> Samples,
        Span<byte> Destination,
        double[] ImpulseResponse,
        double[] StateI,
        double[] StateQ)
    {
        ValidateArguments(ImpulseResponse, StateI, StateQ);

        if ((Samples.Length & 0b11) != 0)
            throw new InvalidOperationException("Длина входного буфера interleaved IQ int16 должна быть кратна 4");

        if (Destination.Length < Samples.Length)
            throw new ArgumentException("Размер буфера назначения должен быть не меньше размера входного буфера", nameof(Destination));

        for (var i = 0; i < Samples.Length; i += 4)
        {
            var i_input = BinaryPrimitives.ReadInt16LittleEndian(Samples.Slice(i, 2));
            var q_input = BinaryPrimitives.ReadInt16LittleEndian(Samples.Slice(i + 2, 2));

            var i_filtered = StateI.FilterSample(ImpulseResponse, i_input);
            var q_filtered = StateQ.FilterSample(ImpulseResponse, q_input);

            BinaryPrimitives.WriteInt16LittleEndian(Destination.Slice(i, 2), ClampToInt16(i_filtered));
            BinaryPrimitives.WriteInt16LittleEndian(Destination.Slice(i + 2, 2), ClampToInt16(q_filtered));
        }
    }

    /// <summary>Отфильтровать raw буфер interleaved IQ int16 little-endian с записью в этот же буфер</summary>
    /// <param name="Samples">Буфер формата I16LE,Q16LE,I16LE,Q16LE</param>
    /// <param name="ImpulseResponse">Массив коэффициентов КИХ-фильтра</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    public static void FilterFIRInterleavedInt16Le(
        this Span<byte> Samples,
        double[] ImpulseResponse,
        double[] StateI,
        double[] StateQ) => FilterFIRInterleavedInt16Le((ReadOnlySpan<byte>)Samples, Samples, ImpulseResponse, StateI, StateQ);

    /// <summary>Отфильтровать raw буфер interleaved IQ int16 little-endian с автоматическим созданием векторов состояния</summary>
    /// <param name="Samples">Входной буфер формата I16LE,Q16LE,I16LE,Q16LE</param>
    /// <param name="ImpulseResponse">Массив коэффициентов КИХ-фильтра</param>
    /// <returns>Отфильтрованный буфер того же формата</returns>
    public static byte[] FilterFIRInterleavedInt16Le(this ReadOnlySpan<byte> Samples, double[] ImpulseResponse)
    {
        ArgumentNullException.ThrowIfNull(ImpulseResponse);

        var result = new byte[Samples.Length];
        Samples.FilterFIRInterleavedInt16Le(result, ImpulseResponse, new double[ImpulseResponse.Length], new double[ImpulseResponse.Length]);
        return result;
    }

    /// <summary>Отфильтровать interleaved IQ int16 LE поток из входного потока в выходной поток</summary>
    /// <param name="Source">Входной поток с данными формата I16LE,Q16LE,I16LE,Q16LE</param>
    /// <param name="Destination">Выходной поток для записи отфильтрованных данных</param>
    /// <param name="Buffer">Рабочий буфер чтения/записи</param>
    /// <param name="ImpulseResponse">Массив коэффициентов КИХ-фильтра</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    /// <returns>Количество записанных байт</returns>
    public static long FilterFIRInterleavedInt16Le(
        this Stream Source,
        Stream Destination,
        byte[] Buffer,
        double[] ImpulseResponse,
        double[] StateI,
        double[] StateQ)
    {
        ArgumentNullException.ThrowIfNull(Source);
        ArgumentNullException.ThrowIfNull(Destination);
        ArgumentNullException.ThrowIfNull(Buffer);
        ValidateArguments(ImpulseResponse, StateI, StateQ);

        if (Buffer.Length < 4)
            throw new InvalidOperationException("Размер рабочего буфера должен быть не меньше 4 байт");

        var frame_buffer = new byte[4];
        var tail = new byte[4];
        var tail_count = 0;
        long written_bytes = 0;

        while (true)
        {
            var readed = Source.Read(Buffer, 0, Buffer.Length);
            if (readed == 0)
                break;

            var offset = 0;
            if (tail_count > 0)
            {
                var need = 4 - tail_count;
                var copy = Math.Min(need, readed);
                Buffer.AsSpan(0, copy).CopyTo(tail.AsSpan(tail_count));
                tail_count += copy;
                offset = copy;

                if (tail_count == 4)
                {
                    ProcessFrame(tail, frame_buffer, ImpulseResponse, StateI, StateQ);
                    Destination.Write(frame_buffer, 0, 4);
                    written_bytes += 4;
                    tail_count = 0;
                }
            }

            var remaining = readed - offset;
            var block_count = remaining & ~0b11;
            if (block_count > 0)
            {
                var block = Buffer.AsSpan(offset, block_count);
                FilterFIRInterleavedInt16Le(block, block, ImpulseResponse, StateI, StateQ);
                Destination.Write(Buffer, offset, block_count);
                written_bytes += block_count;
            }

            var tail_bytes = remaining - block_count;
            if (tail_bytes > 0)
            {
                Buffer.AsSpan(offset + block_count, tail_bytes).CopyTo(tail);
                tail_count = tail_bytes;
            }
        }

        if (tail_count != 0)
            throw new IOException("Общая длина данных во входном потоке interleaved IQ int16 должна быть кратна 4");

        return written_bytes;
    }

    /// <summary>Отфильтровать interleaved IQ int16 LE поток из входного потока в выходной поток с автоматическим созданием векторов состояния</summary>
    /// <param name="Source">Входной поток с данными формата I16LE,Q16LE,I16LE,Q16LE</param>
    /// <param name="Destination">Выходной поток для записи отфильтрованных данных</param>
    /// <param name="Buffer">Рабочий буфер чтения/записи</param>
    /// <param name="ImpulseResponse">Массив коэффициентов КИХ-фильтра</param>
    /// <returns>Количество записанных байт</returns>
    public static long FilterFIRInterleavedInt16Le(this Stream Source, Stream Destination, byte[] Buffer, double[] ImpulseResponse)
    {
        ArgumentNullException.ThrowIfNull(ImpulseResponse);

        return FilterFIRInterleavedInt16Le(Source, Destination, Buffer, ImpulseResponse, new double[ImpulseResponse.Length], new double[ImpulseResponse.Length]);
    }

    /// <summary>Асинхронно отфильтровать interleaved IQ int16 LE поток из входного потока в выходной поток</summary>
    /// <param name="Source">Входной поток с данными формата I16LE,Q16LE,I16LE,Q16LE</param>
    /// <param name="Destination">Выходной поток для записи отфильтрованных данных</param>
    /// <param name="Buffer">Рабочий буфер чтения/записи</param>
    /// <param name="ImpulseResponse">Массив коэффициентов КИХ-фильтра</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    /// <param name="Cancellation">Маркер отмены операции</param>
    /// <returns>Количество записанных байт</returns>
    public static async Task<long> FilterFIRInterleavedInt16LeAsync(
        this Stream Source,
        Stream Destination,
        byte[] Buffer,
        double[] ImpulseResponse,
        double[] StateI,
        double[] StateQ,
        CancellationToken Cancellation = default)
    {
        ArgumentNullException.ThrowIfNull(Source);
        ArgumentNullException.ThrowIfNull(Destination);
        ArgumentNullException.ThrowIfNull(Buffer);
        ValidateArguments(ImpulseResponse, StateI, StateQ);

        if (Buffer.Length < 4)
            throw new InvalidOperationException("Размер рабочего буфера должен быть не меньше 4 байт");

        var frame_buffer = new byte[4];
        var tail = new byte[4];
        var tail_count = 0;
        long written_bytes = 0;

        while (true)
        {
            var readed = await Source.ReadAsync(Buffer, 0, Buffer.Length, Cancellation).ConfigureAwait(false);
            if (readed == 0)
                break;

            var offset = 0;
            if (tail_count > 0)
            {
                var need = 4 - tail_count;
                var copy = Math.Min(need, readed);
                Buffer.AsSpan(0, copy).CopyTo(tail.AsSpan(tail_count));
                tail_count += copy;
                offset = copy;

                if (tail_count == 4)
                {
                    ProcessFrame(tail, frame_buffer, ImpulseResponse, StateI, StateQ);
                    await Destination.WriteAsync(frame_buffer, 0, 4, Cancellation).ConfigureAwait(false);
                    written_bytes += 4;
                    tail_count = 0;
                }
            }

            var remaining = readed - offset;
            var block_count = remaining & ~0b11;
            if (block_count > 0)
            {
                var block = Buffer.AsSpan(offset, block_count);
                FilterFIRInterleavedInt16Le(block, block, ImpulseResponse, StateI, StateQ);
                await Destination.WriteAsync(Buffer, offset, block_count, Cancellation).ConfigureAwait(false);
                written_bytes += block_count;
            }

            var tail_bytes = remaining - block_count;
            if (tail_bytes > 0)
            {
                Buffer.AsSpan(offset + block_count, tail_bytes).CopyTo(tail);
                tail_count = tail_bytes;
            }
        }

        if (tail_count != 0)
            throw new IOException("Общая длина данных во входном потоке interleaved IQ int16 должна быть кратна 4");

        return written_bytes;
    }

    /// <summary>Асинхронно отфильтровать interleaved IQ int16 LE поток из входного потока в выходной поток с автоматическим созданием векторов состояния</summary>
    /// <param name="Source">Входной поток с данными формата I16LE,Q16LE,I16LE,Q16LE</param>
    /// <param name="Destination">Выходной поток для записи отфильтрованных данных</param>
    /// <param name="Buffer">Рабочий буфер чтения/записи</param>
    /// <param name="ImpulseResponse">Массив коэффициентов КИХ-фильтра</param>
    /// <param name="Cancellation">Маркер отмены операции</param>
    /// <returns>Количество записанных байт</returns>
    public static Task<long> FilterFIRInterleavedInt16LeAsync(
        this Stream Source,
        Stream Destination,
        byte[] Buffer,
        double[] ImpulseResponse,
        CancellationToken Cancellation = default)
    {
        ArgumentNullException.ThrowIfNull(ImpulseResponse);

        return FilterFIRInterleavedInt16LeAsync(Source, Destination, Buffer, ImpulseResponse, new double[ImpulseResponse.Length], new double[ImpulseResponse.Length], Cancellation);
    }

    private static void ProcessFrame(ReadOnlySpan<byte> Input, Span<byte> Output, double[] ImpulseResponse, double[] StateI, double[] StateQ)
    {
        var i_input = BinaryPrimitives.ReadInt16LittleEndian(Input[..2]);
        var q_input = BinaryPrimitives.ReadInt16LittleEndian(Input[2..4]);

        var i_filtered = StateI.FilterSample(ImpulseResponse, i_input);
        var q_filtered = StateQ.FilterSample(ImpulseResponse, q_input);

        BinaryPrimitives.WriteInt16LittleEndian(Output[..2], ClampToInt16(i_filtered));
        BinaryPrimitives.WriteInt16LittleEndian(Output[2..4], ClampToInt16(q_filtered));
    }

    private static void ValidateArguments(double[] ImpulseResponse, double[] StateI, double[] StateQ)
    {
        ArgumentNullException.ThrowIfNull(ImpulseResponse);
        ArgumentNullException.ThrowIfNull(StateI);
        ArgumentNullException.ThrowIfNull(StateQ);

        if (ImpulseResponse.Length == 0)
            throw new InvalidOperationException("Размер массива импульсной характеристики должен быть больше 0");

        if (StateI.Length != ImpulseResponse.Length)
            throw new InvalidOperationException("Размер вектора состояния I-компоненты должен совпадать с размером массива импульсной характеристики");

        if (StateQ.Length != ImpulseResponse.Length)
            throw new InvalidOperationException("Размер вектора состояния Q-компоненты должен совпадать с размером массива импульсной характеристики");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static short ClampToInt16(double Value) => Value switch
    {
        > short.MaxValue => short.MaxValue,
        < short.MinValue => short.MinValue,
        _ => (short)Math.Round(Value)
    };
}
