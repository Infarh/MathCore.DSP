using MathCore.DSP.Fourier;

using RND = System.Random;

namespace MathCore.DSP.Signals;

public static partial class Signal
{
    public static class Random
    {
        private static readonly Lazy<RND> __Random = new(() => new());

        public static DigitalSignal Uniform(double dt, int SamplesCount, RND? rnd = null)
        {
            rnd ??= __Random.Value;

            var samples = new double[SamplesCount];
            for (var i = 0; i < SamplesCount; i++)
                samples[i] = rnd.NextDouble() * 2 - 1;

            return new SamplesDigitalSignal(dt, samples);
        }

        public static DigitalSignal Normal(double dt, int SamplesCount, double D = 1, double M = 0, RND? rnd = null)
        {
            rnd ??= __Random.Value;

            var samples = rnd.NextNormal(SamplesCount, D, M);

            return new SamplesDigitalSignal(dt, samples);
        }

        public static DigitalSignal SpectrumBand(double dt, int SamplesCount, (double Fmin, double Fmax) Band, RND? rnd = null) =>
            SpectrumBand(dt, SamplesCount, Band.Fmin, Band.Fmax, rnd);

        public static DigitalSignal SpectrumBand(double dt, int SamplesCount, Interval<double> Band, RND? rnd = null) =>
            SpectrumBand(dt, SamplesCount, Band.Min, Band.Max, rnd);

        public static DigitalSignal SpectrumBand(double dt, int SamplesCount, Interval Band, RND? rnd = null) =>
            SpectrumBand(dt, SamplesCount, Band.Min, Band.Max, rnd);

        public static DigitalSignal SpectrumBand(double dt, int SamplesCount, double Fmin, double Fmax, RND? rnd = null)
        {
            rnd ??= __Random.Value;

            var df = 1 / dt / SamplesCount;

            var spectrum = new Complex[SamplesCount];
            var i_min = (int)(Fmin / df);
            var i_max = Math.Min((int)(Fmax / df), SamplesCount - 1);
            for (var i = i_min; i <= i_max; i++)
            {
                var re = rnd.NextDouble() - 0.5;
                var im = rnd.NextDouble() - 0.5;
                spectrum[i] = new(re, im);
            }

            var samples_complex = spectrum.FastFourierInverse();

            var samples = new double[samples_complex.Length];

            var p = 0d;
            for (var i = 0; i < samples_complex.Length; i++)
            {
                var re = samples_complex[i].Re;
                p += re * re;
                samples[i] = re;
            }

            p = Math.Sqrt(p / samples_complex.Length);
            for (var i = 0; i < samples.Length; i++)
                samples[i] /= p;

            var signal = new SamplesDigitalSignal(dt, samples);

            return signal;
        }
    }
}