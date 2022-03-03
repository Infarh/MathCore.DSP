namespace MathCore.DSP.Signals;

public static partial class Signal
{
    public static class Harmonic
    {
        public static DigitalSignal Sin(double f0, double dt, int SamplesCount) => MathSamplesSignal.Sin(1, f0, 0, dt, SamplesCount);
        public static DigitalSignal Sin(double A, double f0, double dt, int SamplesCount) => MathSamplesSignal.Sin(A, f0, 0, dt, SamplesCount);
        public static DigitalSignal Sin(double A, double f0, double phi0, double dt, int SamplesCount) => MathSamplesSignal.Sin(A, f0, phi0, dt, SamplesCount);

        public static DigitalSignal Cos(double f0, double dt, int SamplesCount) => MathSamplesSignal.Cos(1, f0, 0, dt, SamplesCount);
        public static DigitalSignal Cos(double A, double f0, double dt, int SamplesCount) => MathSamplesSignal.Cos(A, f0, 0, dt, SamplesCount);
        public static DigitalSignal Cos(double A, double f0, double phi0, double dt, int SamplesCount) => MathSamplesSignal.Cos(A, f0, phi0, dt, SamplesCount);
    }
}