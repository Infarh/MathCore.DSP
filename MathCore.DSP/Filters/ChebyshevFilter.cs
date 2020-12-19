﻿using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    public abstract class ChebyshevFilter : AnalogBasedFilter
    {
        /// <inheritdoc />
        protected ChebyshevFilter([NotNull] double[] B, [NotNull] double[] A) : base(B, A) { }
    }
}