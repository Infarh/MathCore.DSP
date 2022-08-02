using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text;
using System.Xml.Linq;

using MathCore;

using static MathCore.Values.CSV;
// ReSharper disable CheckNamespace
namespace System.Diagnostics;

internal static class DebugEx
{

    public static string ThrowIfNullOrEmpty(this string str, [CallerArgumentExpression("str")] string? StrName = null)
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

    public static void ToDebugEnum(this IEnumerable items, [CallerArgumentExpression("items")] string? Name = null)
    {
        if (Name is { Length: > 0 })
            Debug.WriteLine("object[] {0} = {{", Name);
        var i = 0;
        var culture = CultureInfo.InvariantCulture;
        foreach (var item in items)
        {
            if (i > 0)
                Debug.WriteLine(",");

            FormattableString msg = $"            /*[{i,2}]*/ {item}";
            Debug.Write(msg.ToString(culture));

            i++;
        }
        Debug.WriteLine("");
        Debug.WriteLine("}");
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

            FormattableString msg = $"            /*[{i,2}]*/ {item}";
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
        foreach (var z in items)
        {
            if (i > 0)
                Debug.WriteLine(",");

            //var msg = im == 0
            //    ? (FormattableString)
            //      $"            /*[{i,2}]*/  {re:F18}"
            //    : $"            /*[{i,2}]*/ ({re:F18}, {im:F18})";
            var msg = z switch
            {
                (var re, 0) => (FormattableString)$"            /*[{i,2}]*/  {re:F18}",
                (0, var im) => (FormattableString)$"            /*[{i,2}]*/  (0, {im:F18})",
                var (re, im) => (FormattableString)$"            /*[{i,2}]*/ ({re:F18}, {im:F18})"
            };

            Debug.Write(msg.ToString(culture));

            i++;
        }
        Debug.WriteLine("");
        Debug.WriteLine("}");
    }

    public static TestResult ToDebugEnum(this TestResult result, IEnumerable<Complex> items, [CallerArgumentExpression("items")] string? Name = null)
    {
        var log = new StringBuilder(result.LogOutput);

        if (Name is { Length: > 0 })
            log.AppendFormat("Complex[] {0} = {{\r\n", Name);
        var i = 0;
        var culture = CultureInfo.InvariantCulture;
        foreach (var z in items)
        {
            if (i > 0)
                log.AppendLine(",");

            var msg = z switch
            {
                (var re, 0) => (FormattableString)$"            /*[{i,2}]*/  {re:F18}",
                (0, var im) => (FormattableString)$"            /*[{i,2}]*/  (0, {im:F18})",
                var (re, im) => (FormattableString)$"            /*[{i,2}]*/ ({re:F18}, {im:F18})"
            };

            log.Append(msg.ToString(culture));

            i++;
        }
        log.AppendLine("");
        log.AppendLine("}");

        result.LogOutput = log.ToString();
        return result;
    }
}
