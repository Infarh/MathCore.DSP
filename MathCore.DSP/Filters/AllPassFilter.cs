namespace MathCore.DSP.Filters
{
    /// <summary>Всепропускающий фильльт</summary>
    // http://www.dsplib.ru/content/allpass/allpass.html
    public class AllPassFilter : IIR
    {
        public AllPassFilter(double k) : base(new []{ k, 1 }, new [] { 1, k }) { }
    }
}