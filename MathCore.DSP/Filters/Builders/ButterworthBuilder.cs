namespace MathCore.DSP.Filters.Builders;

/// <summary>Построитель фильтра Батерворта</summary>
public readonly ref struct ButterworthBuilder
{
    /// <summary>Тип пропускания полосы частот</summary>
    public FrequencyPassType PassType { get; init; }

    /// <summary>Период дискретизации фильтра</summary>
    public double dt { get; init; }

    /// <summary>Частота дискретизации фильтра</summary>
    public double fd { get => 1 / dt; init => dt = 1 / value; }

    /// <summary>Фильтр нижних частот</summary>
    public ButterworthBuilder LowPass => this with { PassType = FrequencyPassType.LowPass };

    /// <summary>Фильтр верхних частот</summary>
    public ButterworthBuilder HighPass => this with { PassType = FrequencyPassType.HighPass };

    /// <summary>Полосопропускающий фильтр</summary>
    public ButterworthBuilder BandPass => this with { PassType = FrequencyPassType.BandPass };

    /// <summary>Полосозадерживающий фильтр</summary>
    public ButterworthBuilder BandStop => this with { PassType = FrequencyPassType.BandStop };

    /// <summary>Частота пропускания</summary>
    public double? PassFrequency { get; init; }

    /// <summary>Частота заграждения</summary>
    public double? StopFrequency { get; init; }

    /// <summary>Верхняя частота пропускания</summary>
    public double? PassHighFrequency { get; init; }

    /// <summary>Верхняя частота заграждения</summary>
    public double? StopHighFrequency { get; init; }

    /// <summary>Порядок фильтра</summary>
    public int? Order { get; init; }

    /// <summary>Коэффициент пропускания</summary>
    public double? Gp { get; init; }

    /// <summary>Коэффициент заграждения</summary>
    public double? Gs { get; init; }

    /// <summary>Неравномерность в полосе пропускания</summary>
    public double? Rp { get; init; }

    /// <summary>Неравномерность в полосе заграждения</summary>
    public double? Rs { get; init; }
}