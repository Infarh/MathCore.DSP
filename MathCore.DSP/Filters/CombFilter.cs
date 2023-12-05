namespace MathCore.DSP.Filters;

/// <summary>Гребенчатый фильтр</summary>
/// <remarks>Инициализация нового гребенчатого фильтра</remarks>
/// <param name="D">Задержка</param>
public class CombFilter(int D) : FIR(GetImpulseResponse(D))
{
    /// <summary>Импульсная характеристика гребенчатого фильтра</summary>
    /// <param name="DelayLineLength">Длина линии задержки</param>
    /// <returns>Массив значений импульсной характеристики гребенчатого фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Если размер линии задержки меньше 1</exception>
    private static double[] GetImpulseResponse(int DelayLineLength)
    {
        if(DelayLineLength < 1) throw new ArgumentOutOfRangeException(nameof(DelayLineLength), "Длина линии задержки должна быть больше 1");

        var result = new double[DelayLineLength + 1];
        result[0] = 1;
        result[DelayLineLength] = -1;
        return result;
    }
}