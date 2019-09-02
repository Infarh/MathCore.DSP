namespace MathCore.DSP.Filters
{
    public abstract class ButterworthFilter : AnalogBasedFilter
    {
        protected ButterworthFilter(double[] B, double[] A) : base(B, A) { }
    }
}