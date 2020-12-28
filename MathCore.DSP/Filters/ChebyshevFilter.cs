using System.Runtime.Serialization;

using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    /// <summary>Фильтр Чебышева</summary>
    [KnownType(typeof(ChebyshevLowPass))]
    public abstract class ChebyshevFilter : AnalogBasedFilter
    {
        /// <inheritdoc />
        protected ChebyshevFilter([NotNull] double[] B, [NotNull] double[] A) : base(B, A) { }
    }
}