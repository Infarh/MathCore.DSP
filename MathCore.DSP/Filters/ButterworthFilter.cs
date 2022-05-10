using System.Runtime.Serialization;

//using static System.Math;

namespace MathCore.DSP.Filters;

/// <summary>Фильтр Баттерворта</summary>
[KnownType(typeof(ButterworthLowPass))]
[KnownType(typeof(ButterworthBandPass))]
[KnownType(typeof(ButterworthBandStop))]
[KnownType(typeof(ButterworthHighPass))]
public abstract class ButterworthFilter : AnalogBasedFilter
{
    /// <summary>Получить список полюсов нормированного фильтра</summary>
    /// <param name="N">Число полюсов</param>
    /// <param name="EpsP">Затухание фильтра</param>
    /// <param name="W0">Множитель коэффициента затухания</param>
    /// <returns>Массив полюсов нормированного фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Если число полюсов меньше 1</exception>
    protected static IEnumerable<Complex> GetNormPoles(int N, double EpsP, double W0 = 1)
    {
        if (N <= 0) throw new ArgumentOutOfRangeException(nameof(N), N, "Число полюсов должно быть больше 0");

        var r = N % 2; // Нечётность порядка фильтра

        // Радиус окружности размещения полюсов фильтра
        var alpha = W0 * EpsP.Pow(-1d / N);

        // Если порядок фильтра нечётный, то первым добавляем центральный полюс
        if (r != 0) yield return -alpha;
        // Расчёт полюсов
        // Угловой шаг между полюсами
        for (var (i, th0) = (r, Consts.pi05 / N); i < N; i += 2)
        {
            var z = Complex.Exp(alpha, th0 * (i - r + 1));
            yield return z.ComplexConjugate;
            yield return z;
            //var w = th0 * (i - r + 1);
            //var sin = -alpha * Sin(w);
            //var cos = +alpha * Cos(w);
            //yield return new Complex(sin, +cos);
            //yield return new Complex(sin, -cos);
        }
    }

    /// <summary>Инициализация фильтра Баттерворта</summary>
    /// <param name="B">Коэффициенты полинома числителя</param>
    /// <param name="A">Коэффициенты полинома знаменателя</param>
    /// <param name="Spec">Спецификация фильтра</param>
    protected ButterworthFilter(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
}