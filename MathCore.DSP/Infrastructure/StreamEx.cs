using System.Buffers;

namespace System.IO;

internal static class StreamEx
{
    public static int ReadToSpan(this Stream steam, Span<byte> buffer)
    {
        var shared_buffer = ArrayPool<byte>.Shared.Rent(buffer.Length);
        try
        {
            var readed = steam.Read(shared_buffer, 0, buffer.Length);
            if ((uint)readed > buffer.Length)
                throw new IOException($"Длина потока недостаточна для чтения данных в буфер указанной длины {buffer.Length}")
                {
                    Data =
                    {
                        { "buffer.Length", buffer.Length }
                    }
                };

            new Span<byte>(shared_buffer, 0, readed).CopyTo(buffer);
            return readed;
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(shared_buffer);
        }
    }
}