﻿using static System.Math;

namespace MathCore.DSP.Filters;

/// <summary>Цифровой ФНЧ на основе RC-цепочки</summary>
public class ExponentialLowPassRC : IIR
{
    /// <summary>Инициализация нового цифрового ФНЧ на основе RC-цепочки</summary>
    /// <param name="f0">Частота среза</param>
    /// <param name="dt">Период дискретизации</param>
    public ExponentialLowPassRC(double f0, double dt) : this(Exp(-2 * PI * f0 * dt)) { }

    private ExponentialLowPassRC(double z) : base(A: [1, z], B: [0, 1 - z]) { }
}