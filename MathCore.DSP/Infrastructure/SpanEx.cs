using System.Runtime.InteropServices;

namespace MathCore.DSP.Infrastructure;

internal static class SpanEx
{
    public static Span<T> Cast<T>(this Span<byte> buffer) where T : unmanaged => MemoryMarshal.Cast<byte, T>(buffer);
}