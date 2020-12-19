using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    public class ChebyshevLowPass : ChebyshevFilter
    {
        public ChebyshevLowPass([NotNull] double[] B, [NotNull] double[] A) : base(B, A)
        {
        }
    }
}
