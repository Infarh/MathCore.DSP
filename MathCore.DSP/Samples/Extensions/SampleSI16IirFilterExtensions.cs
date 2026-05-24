using System.Runtime.CompilerServices;

namespace MathCore.DSP.Samples.Extensions;

/// <summary>Методы БИХ-фильтрации для последовательностей квадратурных целочисленных отсчётов</summary>
public static class SampleSI16IirFilterExtensions
{
    /// <summary>Отфильтровать поток квадратурных отсчётов БИХ-фильтром</summary>
    /// <param name="Samples">Последовательность входных IQ-отсчётов</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    /// <returns>Последовательность отфильтрованных IQ-отсчётов</returns>
    /// <exception cref="ArgumentNullException">Передана пустая ссылка в одном из параметров</exception>
    /// <exception cref="InvalidOperationException">Размеры коэффициентов или состояний фильтра некорректны</exception>
    public static IEnumerable<SampleSI16> FilterIIR(
        this IEnumerable<SampleSI16> Samples,
        double[] A,
        double[] B,
        double[] StateI,
        double[] StateQ)
    {
        ValidateArguments(Samples, A, B, StateI, StateQ);

        foreach (var sample in Samples)
        {
            var i_filtered = StateI.FilterSample(A, B, sample.I);
            var q_filtered = StateQ.FilterSample(A, B, sample.Q);

            yield return new(ClampToSByte(i_filtered), ClampToSByte(q_filtered));
        }
    }

    /// <summary>Отфильтровать поток квадратурных отсчётов БИХ-фильтром с автоматическим созданием векторов состояния</summary>
    /// <param name="Samples">Последовательность входных IQ-отсчётов</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <returns>Последовательность отфильтрованных IQ-отсчётов</returns>
    /// <exception cref="ArgumentNullException">Передана пустая ссылка в одном из параметров</exception>
    /// <exception cref="InvalidOperationException">Размеры коэффициентов фильтра некорректны</exception>
    public static IEnumerable<SampleSI16> FilterIIR(this IEnumerable<SampleSI16> Samples, double[] A, double[] B)
    {
        ArgumentNullException.ThrowIfNull(A);

        return FilterIIR(Samples, A, B, new double[A.Length], new double[A.Length]);
    }

    /// <summary>Отфильтровать блок квадратурных отсчётов БИХ-фильтром</summary>
    /// <param name="Samples">Блок входных IQ-отсчётов</param>
    /// <param name="Destination">Буфер для отфильтрованных IQ-отсчётов</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    /// <exception cref="ArgumentNullException">Передана пустая ссылка в одном из параметров</exception>
    /// <exception cref="ArgumentException">Размер буфера назначения меньше размера входного блока</exception>
    /// <exception cref="InvalidOperationException">Размеры коэффициентов или состояний фильтра некорректны</exception>
    public static void FilterIIR(
        this ReadOnlySpan<SampleSI16> Samples,
        Span<SampleSI16> Destination,
        double[] A,
        double[] B,
        double[] StateI,
        double[] StateQ)
    {
        ValidateArguments(A, B, StateI, StateQ);

        if (Destination.Length < Samples.Length)
            throw new ArgumentException("Размер буфера назначения должен быть не меньше размера входного блока", nameof(Destination));

        for (var i = 0; i < Samples.Length; i++)
        {
            var sample = Samples[i];

            var i_filtered = StateI.FilterSample(A, B, sample.I);
            var q_filtered = StateQ.FilterSample(A, B, sample.Q);

            Destination[i] = new(ClampToSByte(i_filtered), ClampToSByte(q_filtered));
        }
    }

    /// <summary>Отфильтровать блок квадратурных отсчётов БИХ-фильтром с автоматическим созданием векторов состояния</summary>
    /// <param name="Samples">Блок входных IQ-отсчётов</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <returns>Массив отфильтрованных IQ-отсчётов</returns>
    /// <exception cref="ArgumentNullException">Передана пустая ссылка в одном из параметров</exception>
    /// <exception cref="InvalidOperationException">Размеры коэффициентов фильтра некорректны</exception>
    public static SampleSI16[] FilterIIR(this ReadOnlySpan<SampleSI16> Samples, double[] A, double[] B)
    {
        ArgumentNullException.ThrowIfNull(A);

        var result = new SampleSI16[Samples.Length];
        Samples.FilterIIR(result, A, B, new double[A.Length], new double[A.Length]);
        return result;
    }

    /// <summary>Отфильтровать raw буфер interleaved IQ-отсчётов формата HackRF</summary>
    /// <param name="Samples">Входной буфер формата I,Q,I,Q</param>
    /// <param name="Destination">Выходной буфер формата I,Q,I,Q</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    /// <exception cref="ArgumentException">Размер буфера назначения меньше размера входного буфера</exception>
    /// <exception cref="InvalidOperationException">Длина входного буфера должна быть чётной</exception>
    public static void FilterIIRInterleaved(
        this ReadOnlySpan<byte> Samples,
        Span<byte> Destination,
        double[] A,
        double[] B,
        double[] StateI,
        double[] StateQ)
    {
        ValidateArguments(A, B, StateI, StateQ);

        if ((Samples.Length & 1) != 0)
            throw new InvalidOperationException("Длина входного буфера interleaved IQ должна быть чётной");

        if (Destination.Length < Samples.Length)
            throw new ArgumentException("Размер буфера назначения должен быть не меньше размера входного буфера", nameof(Destination));

        for (var i = 0; i < Samples.Length; i += 2)
        {
            var i_input = unchecked((sbyte)Samples[i]);
            var q_input = unchecked((sbyte)Samples[i + 1]);

            var i_filtered = StateI.FilterSample(A, B, i_input);
            var q_filtered = StateQ.FilterSample(A, B, q_input);

            Destination[i] = unchecked((byte)ClampToSByte(i_filtered));
            Destination[i + 1] = unchecked((byte)ClampToSByte(q_filtered));
        }
    }

    /// <summary>Отфильтровать raw буфер interleaved IQ-отсчётов формата HackRF с записью в этот же буфер</summary>
    /// <param name="Samples">Буфер формата I,Q,I,Q</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    /// <exception cref="InvalidOperationException">Длина входного буфера должна быть чётной</exception>
    public static void FilterIIRInterleaved(
        this Span<byte> Samples,
        double[] A,
        double[] B,
        double[] StateI,
        double[] StateQ) => FilterIIRInterleaved((ReadOnlySpan<byte>)Samples, Samples, A, B, StateI, StateQ);

    /// <summary>Отфильтровать raw буфер interleaved IQ-отсчётов формата HackRF с автоматическим созданием векторов состояния</summary>
    /// <param name="Samples">Входной буфер формата I,Q,I,Q</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <returns>Отфильтрованный буфер формата I,Q,I,Q</returns>
    public static byte[] FilterIIRInterleaved(this ReadOnlySpan<byte> Samples, double[] A, double[] B)
    {
        ArgumentNullException.ThrowIfNull(A);

        var result = new byte[Samples.Length];
        Samples.FilterIIRInterleaved(result, A, B, new double[A.Length], new double[A.Length]);
        return result;
    }

    /// <summary>Отфильтровать interleaved IQ-поток из входного потока в выходной поток</summary>
    /// <param name="Source">Входной поток с данными формата I,Q,I,Q</param>
    /// <param name="Destination">Выходной поток для записи отфильтрованных данных</param>
    /// <param name="Buffer">Рабочий буфер чтения/записи</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    /// <returns>Количество записанных байт</returns>
    /// <exception cref="ArgumentNullException">Передана пустая ссылка в одном из параметров</exception>
    /// <exception cref="InvalidOperationException">Размер буфера меньше 2 байт</exception>
    /// <exception cref="IOException">Общая длина данных во входном потоке нечётная</exception>
    public static long FilterIIRInterleaved(
        this Stream Source,
        Stream Destination,
        byte[] Buffer,
        double[] A,
        double[] B,
        double[] StateI,
        double[] StateQ)
    {
        ArgumentNullException.ThrowIfNull(Source);
        ArgumentNullException.ThrowIfNull(Destination);
        ArgumentNullException.ThrowIfNull(Buffer);
        ValidateArguments(A, B, StateI, StateQ);

        if (Buffer.Length < 2)
            throw new InvalidOperationException("Размер рабочего буфера должен быть не меньше 2 байт");

        var pair_buffer = new byte[2];
        var has_carry = false;
        byte carry_byte = 0;
        long written_bytes = 0;

        while (true)
        {
            var readed = Source.Read(Buffer, 0, Buffer.Length);
            if (readed == 0)
                break;

            var offset = 0;
            if (has_carry)
            {
                var i_filtered = StateI.FilterSample(A, B, unchecked((sbyte)carry_byte));
                var q_filtered = StateQ.FilterSample(A, B, unchecked((sbyte)Buffer[0]));

                pair_buffer[0] = unchecked((byte)ClampToSByte(i_filtered));
                pair_buffer[1] = unchecked((byte)ClampToSByte(q_filtered));

                Destination.Write(pair_buffer, 0, 2);
                written_bytes += 2;

                has_carry = false;
                offset = 1;
            }

            var count = readed - offset;
            var even_count = count & ~1;

            if (even_count > 0)
            {
                var source_span = Buffer.AsSpan(offset, even_count);
                var destination_span = Buffer.AsSpan(offset, even_count);

                FilterIIRInterleaved(source_span, destination_span, A, B, StateI, StateQ);

                Destination.Write(Buffer, offset, even_count);
                written_bytes += even_count;
            }

            if (count != even_count)
            {
                carry_byte = Buffer[offset + even_count];
                has_carry = true;
            }
        }

        if (has_carry)
            throw new IOException("Общая длина данных во входном потоке interleaved IQ должна быть чётной");

        return written_bytes;
    }

    /// <summary>Отфильтровать interleaved IQ-поток из входного потока в выходной поток с автоматическим созданием векторов состояния</summary>
    /// <param name="Source">Входной поток с данными формата I,Q,I,Q</param>
    /// <param name="Destination">Выходной поток для записи отфильтрованных данных</param>
    /// <param name="Buffer">Рабочий буфер чтения/записи</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <returns>Количество записанных байт</returns>
    public static long FilterIIRInterleaved(this Stream Source, Stream Destination, byte[] Buffer, double[] A, double[] B)
    {
        ArgumentNullException.ThrowIfNull(A);

        return FilterIIRInterleaved(Source, Destination, Buffer, A, B, new double[A.Length], new double[A.Length]);
    }

    /// <summary>Асинхронно отфильтровать interleaved IQ-поток из входного потока в выходной поток</summary>
    /// <param name="Source">Входной поток с данными формата I,Q,I,Q</param>
    /// <param name="Destination">Выходной поток для записи отфильтрованных данных</param>
    /// <param name="Buffer">Рабочий буфер чтения/записи</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <param name="StateI">Вектор состояния фильтра синфазной компоненты</param>
    /// <param name="StateQ">Вектор состояния фильтра квадратурной компоненты</param>
    /// <param name="Cancellation">Маркер отмены операции</param>
    /// <returns>Количество записанных байт</returns>
    /// <exception cref="ArgumentNullException">Передана пустая ссылка в одном из параметров</exception>
    /// <exception cref="InvalidOperationException">Размер буфера меньше 2 байт</exception>
    /// <exception cref="IOException">Общая длина данных во входном потоке нечётная</exception>
    public static async Task<long> FilterIIRInterleavedAsync(
        this Stream Source,
        Stream Destination,
        byte[] Buffer,
        double[] A,
        double[] B,
        double[] StateI,
        double[] StateQ,
        CancellationToken Cancellation = default)
    {
        ArgumentNullException.ThrowIfNull(Source);
        ArgumentNullException.ThrowIfNull(Destination);
        ArgumentNullException.ThrowIfNull(Buffer);
        ValidateArguments(A, B, StateI, StateQ);

        if (Buffer.Length < 2)
            throw new InvalidOperationException("Размер рабочего буфера должен быть не меньше 2 байт");

        var pair_buffer = new byte[2];
        var has_carry = false;
        byte carry_byte = 0;
        long written_bytes = 0;

        while (true)
        {
            var readed = await Source.ReadAsync(Buffer, 0, Buffer.Length, Cancellation).ConfigureAwait(false);
            if (readed == 0)
                break;

            var offset = 0;
            if (has_carry)
            {
                var i_filtered = StateI.FilterSample(A, B, unchecked((sbyte)carry_byte));
                var q_filtered = StateQ.FilterSample(A, B, unchecked((sbyte)Buffer[0]));

                pair_buffer[0] = unchecked((byte)ClampToSByte(i_filtered));
                pair_buffer[1] = unchecked((byte)ClampToSByte(q_filtered));

                await Destination.WriteAsync(pair_buffer, 0, 2, Cancellation).ConfigureAwait(false);
                written_bytes += 2;

                has_carry = false;
                offset = 1;
            }

            var count = readed - offset;
            var even_count = count & ~1;

            if (even_count > 0)
            {
                var source_span = Buffer.AsSpan(offset, even_count);
                var destination_span = Buffer.AsSpan(offset, even_count);

                FilterIIRInterleaved(source_span, destination_span, A, B, StateI, StateQ);

                await Destination.WriteAsync(Buffer, offset, even_count, Cancellation).ConfigureAwait(false);
                written_bytes += even_count;
            }

            if (count != even_count)
            {
                carry_byte = Buffer[offset + even_count];
                has_carry = true;
            }
        }

        if (has_carry)
            throw new IOException("Общая длина данных во входном потоке interleaved IQ должна быть чётной");

        return written_bytes;
    }

    /// <summary>Асинхронно отфильтровать interleaved IQ-поток из входного потока в выходной поток с автоматическим созданием векторов состояния</summary>
    /// <param name="Source">Входной поток с данными формата I,Q,I,Q</param>
    /// <param name="Destination">Выходной поток для записи отфильтрованных данных</param>
    /// <param name="Buffer">Рабочий буфер чтения/записи</param>
    /// <param name="A">Массив коэффициентов обратной связи</param>
    /// <param name="B">Массив коэффициентов прямой связи</param>
    /// <param name="Cancellation">Маркер отмены операции</param>
    /// <returns>Количество записанных байт</returns>
    public static Task<long> FilterIIRInterleavedAsync(
        this Stream Source,
        Stream Destination,
        byte[] Buffer,
        double[] A,
        double[] B,
        CancellationToken Cancellation = default)
    {
        ArgumentNullException.ThrowIfNull(A);

        return FilterIIRInterleavedAsync(Source, Destination, Buffer, A, B, new double[A.Length], new double[A.Length], Cancellation);
    }

    private static void ValidateArguments(IEnumerable<SampleSI16> Samples, double[] A, double[] B, double[] StateI, double[] StateQ)
    {
        ArgumentNullException.ThrowIfNull(Samples);
        ValidateArguments(A, B, StateI, StateQ);
    }

    private static void ValidateArguments(double[] A, double[] B, double[] StateI, double[] StateQ)
    {
        ArgumentNullException.ThrowIfNull(A);
        ArgumentNullException.ThrowIfNull(B);
        ArgumentNullException.ThrowIfNull(StateI);
        ArgumentNullException.ThrowIfNull(StateQ);

        if (A.Length == 0)
            throw new InvalidOperationException("Размер массива коэффициентов знаменателя должен быть больше 0");

        if (B.Length == 0)
            throw new InvalidOperationException("Размер массива коэффициентов числителя должен быть больше 0");

        if (A.Length < B.Length)
            throw new InvalidOperationException("Размеры массивов числителя и знаменателя передаточной функции не равны");

        if (StateI.Length != A.Length)
            throw new InvalidOperationException("Размер вектора состояния I-компоненты должен совпадать с размером массива A");

        if (StateQ.Length != A.Length)
            throw new InvalidOperationException("Размер вектора состояния Q-компоненты должен совпадать с размером массива A");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static sbyte ClampToSByte(double Value) => Value switch
    {
        > sbyte.MaxValue => sbyte.MaxValue,
        < sbyte.MinValue => sbyte.MinValue,
        _ => (sbyte)Math.Round(Value)
    };
}