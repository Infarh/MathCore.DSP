using System.Buffers;
using System.Runtime.InteropServices;

namespace MathCore.DSP.Infrastructure;

internal static class MemoryEx
{
    public static Memory<TResult> Cast<TSource, TResult>(Memory<TSource> memory)
        where TSource : unmanaged
        where TResult : unmanaged =>
        typeof(TSource) == typeof(TResult)
            ? (Memory<TResult>)(object)memory
            : new CastMemoryManager<TSource, TResult>(memory).Memory;

    private sealed class CastMemoryManager<TSource, TResult>(Memory<TSource> Source) : MemoryManager<TResult>
        where TSource : unmanaged
        where TResult : unmanaged
    {
        public override Span<TResult> GetSpan() => MemoryMarshal.Cast<TSource, TResult>(Source.Span);

        public override MemoryHandle Pin(int Index = 0) => throw new NotSupportedException();

        public override void Unpin() => throw new NotSupportedException();

        protected override void Dispose(bool disposing) { }
    }

    public static Memory<T> Cast<T>(this Memory<byte> memory) where T : unmanaged => Cast<byte, T>(memory);
}