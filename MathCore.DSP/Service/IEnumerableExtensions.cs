using System.Collections.Generic;

using MathCore;
// ReSharper disable CheckNamespace

namespace System.Linq
{
    public static class IEnumerableExtensions
    {
        public static Complex Multiply(this IEnumerable<Complex> Values) => Values.Aggregate(Complex.Real, (Z, z) => Z * z);
        public static Complex Multiply(this IEnumerable<Complex> Values, Func<Complex, Complex, Complex> Core) => Values.Aggregate(Complex.Real, Core);
        public static Complex Multiply(this IEnumerable<Complex> Values, Func<Complex, Complex> Core) => Values.Aggregate(Complex.Real, (Z, z) => Z * Core(z));
    }
}
