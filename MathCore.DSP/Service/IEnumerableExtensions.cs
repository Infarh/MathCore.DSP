using System.Collections.Generic;

using MathCore;

using static MathCore.Complex;

// ReSharper disable CheckNamespace

namespace System.Linq
{
    // ReSharper disable once InconsistentNaming
    public static class IEnumerableExtensions
    {
        public static Complex Multiply(this IEnumerable<Complex> Values) => Values.Aggregate(Real, (Z, z) => Z * z);
        public static Complex Multiply(this IEnumerable<Complex> Values, Func<Complex, Complex, Complex> Core) => Values.Aggregate(Real, Core);
        public static Complex Multiply(this IEnumerable<Complex> Values, Func<Complex, Complex> Core) => Values.Aggregate(Real, (Z, z) => Z * Core(z));
    }
}
