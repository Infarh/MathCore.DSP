﻿using System.Diagnostics;
using System.Globalization;
using System.Runtime.CompilerServices;

namespace MathCore.DSP.Infrastructure;

internal static class DebugExtensions
{
    public static T ToDebug<T>(this T value, [CallerArgumentExpression("value")] string? Prefix = null)
    {
        if (Prefix is { Length: > 0 })
        {
            FormattableString msg = $"{Prefix} = {value}";
            Debug.WriteLine(msg.ToString(CultureInfo.InvariantCulture));
        }
        else
        {
            FormattableString msg = $"{value}";
            Debug.WriteLine(msg.ToString(CultureInfo.InvariantCulture));
        }
        return value;
    }

    public static void ToDebugEnum<T>(this IEnumerable<T> items, [CallerArgumentExpression("items")] string? Name = null)
    {
        if (Name is { Length: > 0 })
            Debug.WriteLine("{0}[] {1} = {{", typeof(T).Name, Name);
        var i = 0;
        var culture = CultureInfo.InvariantCulture;
        foreach (var item in items)
        {
            if (i > 0)
                Debug.WriteLine(",");

            FormattableString msg = $"/*[{i,2}]*/ {item}";
            Debug.Write(msg.ToString(culture));

            i++;
        }
        Debug.WriteLine("");
        Debug.WriteLine("}");
    }

    public static void ToDebugEnum(this IEnumerable<Complex> items, [CallerArgumentExpression("items")] string? Name = null)
    {
        if (Name is { Length: > 0 })
            Debug.WriteLine("Complex[] {0} = {{", (object)Name);
        var i = 0;
        var culture = CultureInfo.InvariantCulture;
        foreach (var (re, im) in items)
        {
            if (i > 0)
                Debug.WriteLine(",");

            var msg = im == 0
                ? (FormattableString)
                $"            /*[{i,2}]*/  {re:F18}"
                : $"            /*[{i,2}]*/ ({re:F18}, {im:F18})";
            Debug.Write(msg.ToString(culture));

            i++;
        }
        Debug.WriteLine("");
        Debug.WriteLine("}");
    }
}
