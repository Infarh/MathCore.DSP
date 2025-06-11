namespace MathCore.DSP.Filters.Builders;

/// <summary>Строитель эллиптических фильтров</summary>
public readonly ref struct EllipticBuilder
{
    /// <summary>Тип пропускания полосы частот</summary>
    public FrequencyPassType PassType { get; init; }

    /// <summary>Период дискретизации</summary>
    public double dt { get; init; }

    /// <summary>Частота дискретизации</summary>
    public double fd { get => 1 / dt; init => dt = 1 / value; }

    /// <summary>ФНЧ</summary>
    public EllipticBuilder LowPass => this with { PassType = FrequencyPassType.LowPass };

    /// <summary>ФВЧ</summary>
    public EllipticBuilder HighPass => this with { PassType = FrequencyPassType.HighPass };

    /// <summary>ППФ</summary>
    public EllipticBuilder BandPass => this with { PassType = FrequencyPassType.BandPass };

    /// <summary>ПЗФ</summary>
    public EllipticBuilder BandStop => this with { PassType = FrequencyPassType.BandStop };

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

    /// <summary>Коэффициент передачи в полосе пропускания</summary>
    public double? Gp { get; init; }

    /// <summary>Коэффициент передачи в полосе заграждения</summary>
    public double? Gs { get; init; }

    /// <summary>Неравномерность в полосе пропускания (дБ)</summary>
    public double? Rp { get; init; }

    /// <summary>Затухание в полосе заграждения (дБ)</summary>
    public double? Rs { get; init; }
}