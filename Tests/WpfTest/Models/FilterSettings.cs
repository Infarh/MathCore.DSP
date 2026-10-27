#nullable enable
using MathCore.DSP.Filters;

namespace WpfTest.Models;

/// <summary>Параметры синтеза БИХ-фильтра</summary>
public sealed class FilterSettings
{
    /// <summary>Тип полосы фильтра</summary>
    public FilterKind Kind { get; set; } = FilterKind.BandPass;

    /// <summary>Семейство (реализация) фильтра</summary>
    public FilterFamily Family { get; set; } = FilterFamily.Butterworth;

    /// <summary>Частота дискретизации, Гц</summary>
    public double SamplingFrequency { get; set; } = 10_000;

    /// <summary>Нижняя граница полосы пропускания, Гц</summary>
    public double PassLowFrequency { get; set; } = 1_000;

    /// <summary>Верхняя граница полосы пропускания, Гц</summary>
    public double PassHighFrequency { get; set; } = 2_000;

    /// <summary>Нижняя граница полосы заграждения, Гц</summary>
    public double StopLowFrequency { get; set; } = 500;

    /// <summary>Верхняя граница полосы заграждения, Гц</summary>
    public double StopHighFrequency { get; set; } = 2_500;

    /// <summary>Неравномерность АЧХ в полосе пропускания, дБ</summary>
    public double PassRippleDb { get; set; } = 1;

    /// <summary>Подавление в полосе заграждения, дБ</summary>
    public double StopAttenuationDb { get; set; } = 40;

    /// <summary>Частота среза RC-фильтра, Гц</summary>
    public double CutoffFrequency { get; set; } = 1_000;

    /// <summary>Центральная частота RLC-контура, Гц</summary>
    public double CenterFrequency { get; set; } = 1_500;

    /// <summary>Полоса пропускания RLC-контура, Гц</summary>
    public double Bandwidth { get; set; } = 500;

    /// <summary>Частота Найквиста, Гц</summary>
    public double NyquistFrequency => SamplingFrequency / 2;

    /// <summary>Копия параметров</summary>
    /// <returns>Новый экземпляр с теми же значениями</returns>
    public FilterSettings Clone() => (FilterSettings)MemberwiseClone();
}

/// <summary>Коэффициенты полиномов БИХ-фильтра</summary>
/// <param name="SamplingFrequency">Частота дискретизации, Гц</param>
/// <param name="B">Коэффициенты полинома числителя</param>
/// <param name="A">Коэффициенты полинома знаменателя</param>
/// <param name="Description">Описание фильтра</param>
public sealed record FilterCoefficients(double SamplingFrequency, double[] B, double[] A, string? Description = null)
{
    /// <summary>Порядок фильтра</summary>
    public int Order => Math.Max(A.Length, B.Length) - 1;
}

/// <summary>Результат синтеза БИХ-фильтра</summary>
/// <param name="Filter">Синтезированный фильтр</param>
/// <param name="Settings">Параметры синтеза</param>
/// <param name="Description">Описание фильтра</param>
public sealed record FilterDesignResult(IIR Filter, FilterSettings Settings, string Description)
{
    /// <summary>Порядок фильтра</summary>
    public int Order => Filter.Order;

    /// <summary>Коэффициенты полиномов фильтра</summary>
    /// <returns>Коэффициенты полинома числителя и знаменателя</returns>
    public FilterCoefficients GetCoefficients() =>
        new(Settings.SamplingFrequency, Filter.B.ToArray(), Filter.A.ToArray(), Description);
}

/// <summary>Исключение при некорректных параметрах синтеза фильтра</summary>
public sealed class FilterDesignException(string Message) : Exception(Message);
