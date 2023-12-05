namespace MathCore.DSP.Filters;

/// <summary>Параллельный комбинационный фильтр</summary>
/// <remarks>Инициализация нового параллельного комбинационного фильтра</remarks>
/// <param name="Filter1">Первый фильтр в комбинации</param>
/// <param name="Filter2">Второй фильтр в комбинации</param>
public class ParallelFilter(Filter Filter1, Filter Filter2) : CombinationFilter(Filter1, Filter2)
{
    public override double Process(double Sample) => Filter1.Process(Sample / 2) + Filter2.Process(Sample / 2);

    public override void Reset()
    {
        Filter1.Reset();
        Filter2.Reset();
    }

    public override Complex FrequencyResponse(double f) => (Filter1.FrequencyResponse(f) + Filter2.FrequencyResponse(f)) / 2;
}