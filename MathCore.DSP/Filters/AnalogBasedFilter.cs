using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    /// <summary>Цифровой фильтр на основе аналогового прототипа</summary>
    public abstract class AnalogBasedFilter : IIR
    {
        /// <summary>Инициализация параметров цифрового фильтра на базе аналогового прототипа</summary>
        /// <param name="B">Коэффициенты полинома числителя</param>
        /// <param name="A">Коэффициенты полинома знаменателя</param>
        protected AnalogBasedFilter([NotNull] double[] B, [NotNull] double[] A) : base(B, A) { }
    }
}