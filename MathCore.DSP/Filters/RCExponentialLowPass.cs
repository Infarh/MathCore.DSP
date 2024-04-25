using static System.Math;

namespace MathCore.DSP.Filters;

/// <summary>Цифровой ФНЧ на основе RC-цепочки</summary>
public class RCExponentialLowPass : IIR
{
    /// <summary>Инициализация нового цифрового ФНЧ на основе RC-цепочки</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="f0">Частота среза</param>
    public RCExponentialLowPass(double dt, double f0) : this(Exp(-2 * PI * f0 * dt)) { }

    private RCExponentialLowPass(double z) : base(A: [1, z], B: [0, 1 - z]) { }
}