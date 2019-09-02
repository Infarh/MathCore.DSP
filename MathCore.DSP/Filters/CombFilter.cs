using System;

namespace MathCore.DSP.Filters
{
    /// <summary>Гребенчатый фильтр</summary>
    public class CombFilter : FIR
    {
        private static double[] GetImpulseResponse(int D)
        {
            if(D < 1) throw new ArgumentOutOfRangeException(nameof(D), "Длина линии задержки должна быть больше 1");

            var result = new double[D + 1];
            result[0] = 1;
            result[D] = -1;
            return result;
        }

        /// <summary>Новый гребенчатый фильтр</summary>
        /// <param name="D">Задержка</param>
        public CombFilter(int D) : base(GetImpulseResponse(D)) { }
    }
}