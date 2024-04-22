using static System.Math;

namespace MathCore.DSP.Filters;

/// <summary>Полосозадерживающий RLC-фильтр</summary>
public class RLCBandStop : IIR
{
    /// <summary>Инициализация нового экземпляра <see cref="RLCBandStop"/></summary>
    /// <param name="f0">Частота резонанса</param>
    /// <param name="DeltaF">Полоса частот по уровню 0.707</param>
    /// <param name="dt">Период дискретизации</param>
    public RLCBandStop(double f0, double DeltaF, double dt) : this(Tan(PI * f0 * dt), PI * DeltaF * dt) { }

    /// <summary>Инициализация нового экземпляра <see cref="RLCBandStop"/></summary>
    /// <param name="w0">Обобщённая центральная частота</param>
    /// <param name="dw">Обобщённая полоса частот</param>
    private RLCBandStop(double w0, double dw)
        : base(
            A: [w0 * w0 + dw + 1, 2 * (w0 * w0 - 1), w0 * w0 - dw + 1],
            B: [w0*w0 + 1, 2 * (w0 * w0 - 1), w0 * w0 + 1])
    { }
}