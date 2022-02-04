using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace MathCore.DSP.Infrastructure;

internal static class MemoryEx
{
    public static Memory<TResult> Cast<TSource, TResult>(Memory<TSource> memory)
        where TSource : unmanaged
        where TResult : unmanaged
    {
        if (typeof(TSource) == typeof(TResult)) 
            return (Memory<TResult>)(object)memory;

        return new CastMemoryManager<TSource, TResult>(memory).Memory;
    }

    private sealed class CastMemoryManager<TSource, TResult> : MemoryManager<TResult>
        where TSource : unmanaged
        where TResult : unmanaged
    {
        private readonly Memory<TSource> _Source;

        public CastMemoryManager(Memory<TSource> Source) => _Source = Source;

        public override Span<TResult> GetSpan() => MemoryMarshal.Cast<TSource, TResult>(_Source.Span);

        public override MemoryHandle Pin(int Index = 0) => throw new NotSupportedException();

        public override void Unpin() => throw new NotSupportedException();

        protected override void Dispose(bool disposing) { }
    }

    public static Memory<T> Cast<T>(this Memory<byte> memory) where T : unmanaged
    {
        return Cast<byte, T>(memory);
    }
}