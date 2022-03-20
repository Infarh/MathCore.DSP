using System.Runtime.Serialization;

using static System.Math;

namespace MathCore.DSP.Filters;

/// <summary>Фильтр Баттерворта</summary>
[KnownType(typeof(ButterworthLowPass))]
[KnownType(typeof(ButterworthBandPass))]
[KnownType(typeof(ButterworthHighPass))]
public abstract class ButterworthFilter : AnalogBasedFilter
{
    /// <summary>Получить список полюсов нормированного фильтра</summary>
    /// <param name="N">Число полюсов</param>
    /// <param name="EpsP">Затухание фильтра</param>
    /// <returns>Массив полюсов нормированного фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Если число полюсов меньше 1</exception>
    protected static IEnumerable<Complex> GetNormPoles(int N, double EpsP)
    {
        if (N <= 0) throw new ArgumentOutOfRangeException(nameof(N), N, "Число полюсов должно быть больше 0");

        var r = N % 2; // Нечётность порядка фильтра

        // Радиус окружности размещения полюсов фильтра
        var alpha = EpsP.Pow(-1d / N);

        // Угловой шаг между полюсами
        var dth = PI / N;

        // Если порядок фильтра нечётный, то первым добавляем центральный полюс
        if (r != 0) yield return -alpha;
        // Расчёт полюсов
        for (var i = r; i < N; i += 2)
        {
            var w = dth * (i + 1 - r - 0.5);
            var sin = -alpha * Sin(w);
            var cos = alpha * Cos(w);
            yield return new Complex(sin, cos);
            yield return new Complex(sin, -cos);
        }
    }

    /// <summary>Инициализация фильтра Баттерворта</summary>
    /// <param name="B">Коэффициенты полинома числителя</param>
    /// <param name="A">Коэффициенты полинома знаменателя</param>
    /// <param name="Spec">Спецификация фильтра</param>
    protected ButterworthFilter(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
}