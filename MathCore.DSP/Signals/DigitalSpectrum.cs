namespace MathCore.DSP.Signals
{
    public class DigitalSpectrum
    {
        private double _df;
        private double _f0;

        public double fd => _df;
        public double f0 => _f0;

        protected DigitalSpectrum(double df, double f0) { }
    }
}
