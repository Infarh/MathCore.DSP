using System.Runtime.Serialization;

using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    /// <summary>Фильтр Баттерворта</summary>
    [KnownType(typeof(ButterworthLowPass))]
    public abstract class ButterworthFilter : AnalogBasedFilter
    {
        /// <summary>Инициализация фильтра Баттерворта</summary>
        /// <param name="B">Коэффициенты полинома числителя</param>
        /// <param name="A">Коэффициенты полинома знаменателя</param>
        protected ButterworthFilter([NotNull] double[] B, [NotNull] double[] A) : base(B, A) { }
    }
}