namespace MathCore.DSP.Filters;

/// <summary>Последовательный комбинационный фильтр</summary>
/// <remarks>Инициализация нового последовательного комбинационного фильтра</remarks>
/// <param name="Filter1">Первый фильтр в комбинации</param>
/// <param name="Filter2">Второй фильтр в комбинации</param>
public class SerialFilter(Filter Filter1, Filter Filter2) : CombinationFilter(Filter1, Filter2)
{
    public override double Process(double Sample) => Filter2.Process(Filter1.Process(Sample));

    public override void Reset()
    {
        Filter1.Reset();
        Filter2.Reset();
    }

    public override Complex FrequencyResponse(double f) => Filter1.FrequencyResponse(f) * Filter2.FrequencyResponse(f);
}