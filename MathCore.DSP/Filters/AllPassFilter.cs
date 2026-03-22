// ReSharper disable ArgumentsStyleOther
namespace MathCore.DSP.Filters;

/// <summary>Всепропускающий фильтр</summary>
/// <remarks>Инициализация нового экземпляра <see cref="AllPassFilter"/></remarks>
/// <param name="k">Коэффициент передачи</param>
// http://www.dsplib.ru/content/allpass/allpass.html
public class AllPassFilter(double k) : IIR(B: [k, 1], A: [1, k]);