using System;
using System.Runtime.Serialization;

using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    /// <summary>Фильтр Баттерворта</summary>
    [KnownType(typeof(ButterworthLowPass))]
    public abstract class ButterworthFilter : AnalogBasedFilter
    {
        protected static Complex[] GetNormPoles(int N, double EpsP)
        {
            var r = N % 2; // Нечётность порядка фильтра

            // Радиус окружности размещения полюсов фильтра
            var alpha = EpsP.Pow(-1d / N);

            // Угловой шаг между полюсами
            var dth = Math.PI / N;

            // Массив полюсов фильтра
            var poles = new Complex[N];
            // Если порядок фильтра нечётный, то первым добавляем центральный полюс
            if (r != 0) poles[0] = -alpha;
            // Расчёт полюсов
            for (var i = r; i < poles.Length; i += 2)
            {
                var w = dth * (i + 1 - r - 0.5);
                var sin = -alpha * Math.Sin(w);
                var cos = alpha * Math.Cos(w);
                poles[i] = new Complex(sin, cos);
                poles[i + 1] = new Complex(sin, -cos);
            }

            return poles;
        }

        /// <summary>Инициализация фильтра Баттерворта</summary>
        /// <param name="B">Коэффициенты полинома числителя</param>
        /// <param name="A">Коэффициенты полинома знаменателя</param>
        protected ButterworthFilter([NotNull] double[] B, [NotNull] double[] A) : base(B, A) { }
    }
}