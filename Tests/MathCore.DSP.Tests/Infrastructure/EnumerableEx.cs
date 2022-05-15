using System.Collections.Generic;

namespace System.Linq;

internal static class EnumerableEx
{
    public static IEnumerable<T[]> SelectGroup<T>(this IEnumerable<T> items, int GroupSize)
    {
        var buffer = new T[GroupSize];

        var index = 0;
        foreach (var item in items)
        {
            buffer[index] = item;
            index++;
            if (index != GroupSize) continue;

            index = 0;
            yield return buffer;
            buffer = new T[GroupSize];
        }

        if (index == 0) 
            yield break;

        System.Array.Resize(ref buffer, index);
        yield return buffer;
    }
}