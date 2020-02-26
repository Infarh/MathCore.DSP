using System;
using System.Collections.Generic;
using MathCore.Annotations;

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

        [NotNull]
        public static EnumerableSignal Sin(double A, double f0, double phi0, double dt, int SamplesCount = -1)
        {
            var w0 = 2 * Math.PI * f0;
            return new EnumerableSignal(dt, GetSamples(t => A * Math.Sin(w0 * t + phi0), dt, SamplesCount));
        }

        [NotNull]
        public static EnumerableSignal Cos(double A, double f0, double phi0, double dt, int SamplesCount)
        {
            var w0 = 2 * Math.PI * f0;
            return new EnumerableSignal(dt, GetSamples(t => A * Math.Sin(w0 * t + phi0), dt, SamplesCount));
        }

        
    }
}