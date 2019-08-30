using System;

namespace MathCore.DSP.Filters
{
    /// <summary>Цифровой ФНЧ на основе R-C цепочки</summary>
    public class ExponetialLowPassRC : IIR
    {
        /// <summary>Инициализация нового цифрового ФНЧ на основе RC-цепочики</summary>
        /// <param name="f0">Частота среза</param>
        /// <param name="dt">Период дискретизации</param>
        public ExponetialLowPassRC(double f0, double dt) : this(Math.Exp(-2 * Math.PI * f0 * dt)) { }
        private ExponetialLowPassRC(double z) : base(A: new[] { 1, z }, B: new[] { 0, 1 - z }) { }
    }
}
