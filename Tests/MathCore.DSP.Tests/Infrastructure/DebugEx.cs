using System.Collections.Generic;
using System.Xml.Linq;

namespace System.Diagnostics;

internal static class DebugEx
{
    public static T ToDebug<T>(this T value)
    {
        Debug.WriteLine(value);
        return value;
    }

    public static T ToDebug<T>(this T value, string Prefix)
    {
        Debug.WriteLine("{0}:{1}", Prefix, value);
        return value;
    }

    public static void ToDebugEnumerable<T>(this IEnumerable<T> items, string? Name = null)
    {
        Debug.WriteLine("{0}:", (object)(Name ?? "Items"));
        var i = 0;
        foreach (var item in items)
            Debug.WriteLine("[{0,4}] {1}", i++, item);
    }
}
