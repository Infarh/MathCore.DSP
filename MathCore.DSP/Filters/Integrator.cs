namespace MathCore.DSP.Filters
{
    // Url: http://www.dsplib.ru/content/cic/cic.html
    public class Integrator : IIR
    {
        public Integrator() : base(new[] { 1d }, new[] { 1d, -1d }) { }
    }
}