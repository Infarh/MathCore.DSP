using System.Globalization;

namespace MathCore;

internal static class ComplexExtensions
{
    public static string ToTuplesString(this IEnumerable<Complex> values) => values
       .Select(v => (FormattableString)$"({v.Re:F18}, {v.Im:F18})")
       .Select(s => s.ToString(CultureInfo.InvariantCulture))
       .ToSeparatedStr(",\r\n");
}
