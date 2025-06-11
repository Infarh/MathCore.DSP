namespace MathCore.DSP.Filters.Builders;

/// <summary>Тип пропускания полосы частот</summary>
public enum FrequencyPassType
{
    /// <summary>ФНЧ</summary>
    LowPass,
    /// <summary>ФВЧ</summary>
    HighPass,
    /// <summary>ППФ</summary>
    BandPass,
    /// <summary>ПЗФ</summary>
    BandStop
}