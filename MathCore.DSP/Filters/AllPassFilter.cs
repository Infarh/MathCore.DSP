// ReSharper disable ArgumentsStyleOther
namespace MathCore.DSP.Filters
{
    /// <summary>Всепропускающий фильтр</summary>
    // http://www.dsplib.ru/content/allpass/allpass.html
    public class AllPassFilter : IIR
    {
        /// <summary>Инициализация нового экземпляра <see cref="AllPassFilter"/></summary>
        /// <param name="k">Коэффициент передачи</param>
        public AllPassFilter(double k) : base(B: new []{ k, 1 }, A: new [] { 1, k }) { }
    }
}