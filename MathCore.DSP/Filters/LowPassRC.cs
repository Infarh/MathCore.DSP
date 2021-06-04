using static System.Math;

namespace MathCore.DSP.Filters
{
    /// <summary>Цифровой ФНЧ на основе RC-цепочки с билинейным преобразованием</summary>
    public class LowPassRC : IIR
    {
        /// <summary>Инициализация нового цифрового ФНЧ на основе RC-цепочки с билинейным преобразованием</summary>
        /// <param name="f0">Частота среза</param>
        /// <param name="dt">Период дискретизации</param>
        public LowPassRC(double f0, double dt) : this(1 / Tan(PI * f0 * dt)) { } //Math.PI * f0 * dt
        //public LowPassRC(double f0, double dt) : this(Math.PI * f0 * dt) { }

        //private LowPassRC(double w) : base(A: new[] { 1, (1 - w) / (1 + w) }, B: new[] { w / (1 + w), w / (1 + w) }) { }
        private LowPassRC(double w) : base(A: new[] { 1 + w, 1 - w }, B: new[] { 1d, 1d }) { }
    }
}