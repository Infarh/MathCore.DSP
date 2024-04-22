using static System.Math;

namespace MathCore.DSP.Filters;

/// <summary>Цифровой ФНЧ на основе RC-цепочки с билинейным преобразованием</summary>
public class RCLowPass : IIR
{
    /// <summary>Инициализация нового цифрового ФНЧ на основе RC-цепочки с билинейным преобразованием</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="f0">Частота среза</param>
    public RCLowPass(double dt, double f0) : this(1 / Tan(PI * f0 * dt)) { }

    //private RCLowPass(double w) : base(A: [1, (1 - w) / (1 + w)], B: [w / (1 + w), w / (1 + w)]) { }
    private RCLowPass(double w) : base(A: [1 + w, 1 - w], B: [1d, 1d]) { }
}