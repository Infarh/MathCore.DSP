namespace MathCore.DSP.Filters;

/// <summary>Комбинационный фильтр</summary>
public abstract class CombinationFilter : Filter
{
    /// <summary>Первый фильтр в комбинации</summary>
    public Filter Filter1 { get; }

    /// <summary>Второй фильтр в комбинации</summary>
    public Filter Filter2 { get; }

    /// <summary>Инициализация нового комбинационного фильтра</summary>
    /// <param name="Filter1">Первый фильтр в комбинации</param>
    /// <param name="Filter2">Второй фильтр в комбинации</param>
    protected CombinationFilter(Filter Filter1, Filter Filter2)
    {
        this.Filter1 = Filter1.NotNull();
        this.Filter2 = Filter2.NotNull();
    }
}