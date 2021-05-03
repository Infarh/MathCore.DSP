using System;

using MathCore.Annotations;

namespace MathCore.DSP.Signals
{
    public static class MathSamplesSignal
    {
        [NotNull]
        public static SamplesDigitalSignal Sin(double A, double f0, double phi0, double dt, int SamplesCount)
        {
            var samples = new double[SamplesCount];
            var w0 = 2 * Math.PI * f0 * dt;
            for (var i = 0; i < SamplesCount; i++)
                samples[i] = A * Math.Sin(w0 * i + phi0);
            return new SamplesDigitalSignal(dt, samples);
        }

        [NotNull]
        public static SamplesDigitalSignal Cos(double A, double f0, double phi0, double dt, int SamplesCount)
        {
            var samples = new double[SamplesCount];
            var w0 = 2 * Math.PI * f0 * dt;
            for (var i = 0; i < SamplesCount; i++)
                samples[i] = A * Math.Cos(w0 * i + phi0);
            return new SamplesDigitalSignal(dt, samples);
        }
    }
}
