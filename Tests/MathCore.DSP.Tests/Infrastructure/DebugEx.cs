using System.Collections;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text;

using MathCore;

// ReSharper disable CheckNamespace
namespace System.Diagnostics;

internal static class DebugEx
{

    public static string ThrowIfNullOrEmpty(this string str, [CallerArgumentExpression(nameof(str))] string? StrName = null)
    {
        if (string.IsNullOrEmpty(str))
            throw new ArgumentNullException(StrName, "String is null or empty");
        return str;
    }

    //public static T ToDebug<T>(this T value)
    //{
    //    Debug.WriteLine(value);
    //    return value;
    //}

    public static T ToDebug<T>(this T value, [CallerArgumentExpression(nameof(value))] string? Prefix = null)
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

    private static int Log10Int(int x) => (int)Math.Log10(x);

    public static void ToDebugEnum(this IEnumerable items, [CallerArgumentExpression(nameof(items))] string? Name = null)
    {
        string? pad_str = null;
        if (Name is { Length: > 0 })
        {
            Debug.WriteLine("object[] {0} =", (object)Name);
            Debug.WriteLine("[");
            pad_str = "    ";
        }
        var m = items is ICollection { Count: var items_count }
            ? Log10Int(items_count) + 1
            : 2;
        var i = 0;
        var culture = CultureInfo.InvariantCulture;
        foreach (var item in items)
        {
            if (i > 0)
                Debug.WriteLine(",");

            FormattableString msg = $"{pad_str}/*[{i.ToString().PadLeft(m)}]*/ {item}";
            Debug.Write(msg.ToString(culture));

            i++;
        }

        if (pad_str is null) return;

        Debug.WriteLine("");
        Debug.WriteLine("]");
    }

    public static void ToDebugEnum<T>(this IEnumerable<T> items, [CallerArgumentExpression(nameof(items))] string? Name = null)
    {
        string? pad_str = null;
        if (Name is { Length: > 0 })
        {
            Debug.WriteLine("{0}[] {1} =", typeof(T).Name, Name);
            Debug.WriteLine("[");
            pad_str = "    ";
        }
        var m = items is ICollection { Count: var items_count }
            ? Log10Int(items_count) + 1
            : 2;
        var i = 0;
        var culture = CultureInfo.InvariantCulture;
        foreach (var item in items)
        {
            if (i > 0)
                Debug.WriteLine(",");

            FormattableString msg = $"{pad_str}/*[{i.ToString().PadLeft(m)}]*/ {item}";
            Debug.Write(msg.ToString(culture));

            i++;
        }

        if (pad_str is null) return;

        Debug.WriteLine("");
        Debug.WriteLine("]");
    }

    public static void ToDebugEnum(this IEnumerable<Complex> items, [CallerArgumentExpression(nameof(items))] string? Name = null)
    {
        string? pad_str = null;
        if (Name is { Length: > 0 })
        {
            Debug.WriteLine("Complex[] {0} =", (object)Name);
            Debug.WriteLine("[");
            pad_str = "    ";
        }
        var m = items is ICollection { Count: var items_count }
            ? Log10Int(items_count) + 1
            : 2;
        var i = 0;
        var culture = CultureInfo.InvariantCulture;
        foreach (var z in items)
        {
            if (i > 0)
                Debug.WriteLine(",");

            var msg = z switch
            {
                (var re, 0) => (FormattableString)$"{pad_str}/*[{i.ToString().PadLeft(m)}]*/  {re:F18}",
                (0, var im) => (FormattableString)$"{pad_str}/*[{i.ToString().PadLeft(m)}]*/  (0, {im:F18})",
                var (re, im) => (FormattableString)$"{pad_str}/*[{i.ToString().PadLeft(m)}]*/ ({re:F18}, {im:F18})"
            };

            Debug.Write(msg.ToString(culture));

            i++;
        }

        if (pad_str is null) return;

        Debug.WriteLine("");
        Debug.WriteLine("]");
    }

    public static TestResult ToDebugEnum(this TestResult result, IEnumerable<Complex> items, [CallerArgumentExpression(nameof(items))] string? Name = null)
    {
        var log = new StringBuilder(result.LogOutput);

        string? pad_str = null;
        if (Name is { Length: > 0 })
        {
            log.AppendFormat("Complex[] ").Append(Name).AppendLine(" =");
            log.Append('[').AppendLine();
            pad_str = "    ";
        }
        var m = items is ICollection { Count: var items_count }
            ? Log10Int(items_count) + 1
            : 2;
        var i = 0;
        var culture = CultureInfo.InvariantCulture;
        foreach (var z in items)
        {
            if (i > 0)
                log.Append(',').AppendLine();

            var msg = z switch
            {
                (var re, 0) => (FormattableString)$"{pad_str}/*[{i.ToString().PadLeft(m)}]*/  {re:F18}",
                (0, var im) => (FormattableString)$"{pad_str}/*[{i.ToString().PadLeft(m)}]*/  (0, {im:F18})",
                var (re, im) => (FormattableString)$"{pad_str}/*[{i.ToString().PadLeft(m)}]*/ ({re:F18}, {im:F18})"
            };

            log.Append(msg.ToString(culture));

            i++;
        }

        if (pad_str is not null)
        {
            log.AppendLine("");
            log.AppendLine("]");
        }

        result.LogOutput = log.ToString();
        return result;
    }
}
