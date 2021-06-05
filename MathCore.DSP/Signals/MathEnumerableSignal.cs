using System;
using System.Collections.Generic;

namespace MathCore.DSP.Signals
{
    public static class MathEnumerableSignal
    {
        private static IEnumerable<double> GetSamples(Func<double, double> f, double dt, int SamplesCount)
        {
            var x = 0d;
            var infinity = SamplesCount < 0;
            while (infinity || SamplesCount-- > 0)
                yield return f(x += dt);
        }

        public static EnumerableSignal Sin(double A, double f0, double phi0, double dt, int SamplesCount = -1)
        {
            var w0 = 2 * Math.PI * f0;
            return new EnumerableSignal(dt, GetSamples(t => A * Math.Sin(w0 * t + phi0), dt, SamplesCount));
        }

        public static EnumerableSignal Cos(double A, double f0, double phi0, double dt, int SamplesCount)
        {
            var w0 = 2 * Math.PI * f0;
            return new EnumerableSignal(dt, GetSamples(t => A * Math.Sin(w0 * t + phi0), dt, SamplesCount));
        }

        public static EnumerableSignal Const(double A, double dt, int SamplesCount) => new(dt, GetSamples(_ => A, dt, SamplesCount));
    }
}