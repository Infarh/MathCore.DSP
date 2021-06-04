using System;

namespace MathCore.DSP.Filters
{
    /// <summary>Гребенчатый фильтр</summary>
    public class CombFilter : FIR
    {
        /// <summary>Импульсная характеристика гребенчатого фильтра</summary>
        /// <param name="DelayLineLength">Длина линии задержки</param>
        /// <returns>Массив значений импульсной характеристики гребенчатого фильтра</returns>
        private static double[] GetImpulseResponse(int DelayLineLength)
        {
            if(DelayLineLength < 1) throw new ArgumentOutOfRangeException(nameof(DelayLineLength), "Длина линии задержки должна быть больше 1");

            var result = new double[DelayLineLength + 1];
            result[0] = 1;
            result[DelayLineLength] = -1;
            return result;
        }

        /// <summary>Инициализация нового гребенчатого фильтра</summary>
        /// <param name="D">Задержка</param>
        public CombFilter(int D) : base(GetImpulseResponse(D)) { }
    }
}