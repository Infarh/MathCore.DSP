namespace MathCore.DSP.Filters
{
    public abstract class AnalogBasedFilter : IIR
    {
        protected AnalogBasedFilter(double[] B, double[] A) : base(B, A) { }
    }
}