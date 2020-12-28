using System.Runtime.Serialization;

#nullable enable
namespace MathCore.DSP.Filters
{
    /// <summary>Эллиптический фильтр</summary>
    [KnownType(typeof(EllipticLowPass))]
    public abstract class EllipticFilter : AnalogBasedFilter
    {
        /// <summary>Инициализация нового эллиптического фильтра</summary>
        /// <param name="B">Коэффициенты полинома числителя</param>
        /// <param name="A">Коэффициенты полинома знаменателя</param>
        protected EllipticFilter(double[] B, double[] A) : base(B, A) { }
    }
}
