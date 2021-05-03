using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    public class EllipticBandPass : EllipticFilter
    {
        private EllipticBandPass([NotNull] double[] B, [NotNull] double[] A) : base(B, A) { }
    }
}