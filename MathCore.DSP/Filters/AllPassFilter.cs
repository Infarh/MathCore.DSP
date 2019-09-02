namespace MathCore.DSP.Filters
{
    /// <summary>Всепропускающий фильльт</summary>
    // http://www.dsplib.ru/content/allpass/allpass.html
    public class AllPassFilter : IIR
    {
        /// <summary>Инициализация нового экземпляра <see cref="AllPassFilter"/></summary>
        /// <param name="k">Коэффициент передачи</param>
        public AllPassFilter(double k) : base(new []{ k, 1 }, new [] { 1, k }) { }
    }
}