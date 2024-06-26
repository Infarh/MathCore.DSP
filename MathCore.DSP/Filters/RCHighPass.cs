﻿using static System.Math;

namespace MathCore.DSP.Filters;

/// <summary>Цифровой ФВЧ на основе RС-цепочки с билинейным преобразованием</summary>
public class RCHighPass : IIR
{
    /// <summary>Инициализация нового цифрового ФВЧ на основе RС-цепочки с билинейным преобразованием</summary>
    /// <param name="f0">Частота среза</param>
    /// <param name="dt">Период дискретизации</param>
    public RCHighPass(double f0, double dt) : this(Tan(PI * f0 * dt)) { }
    //private RCHighPass(double q) : base(A: new[] { 1, (1 - q) / (1 + q) }, B: new[] { 1 / (1 + q), -1 / (1 + q) }) { }
    private RCHighPass(double w) : base(A: [w + 1, w - 1], B: [1d, -1d]) { }
}