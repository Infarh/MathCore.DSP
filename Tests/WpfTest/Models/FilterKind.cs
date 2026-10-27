#nullable enable
namespace WpfTest.Models;

/// <summary>Тип полосы пропускания БИХ-фильтра</summary>
public enum FilterKind
{
    /// <summary>Фильтр нижних частот</summary>
    LowPass,
    /// <summary>Фильтр верхних частот</summary>
    HighPass,
    /// <summary>Полосовой фильтр</summary>
    BandPass,
    /// <summary>Полосозаграждающий фильтр</summary>
    BandStop,
}

/// <summary>Описание типа фильтра для пользовательского интерфейса</summary>
/// <param name="Kind">Тип фильтра</param>
/// <param name="Title">Полное название</param>
/// <param name="ShortTitle">Короткое обозначение</param>
/// <param name="Description">Пояснение</param>
public sealed record FilterKindInfo(FilterKind Kind, string Title, string ShortTitle, string Description)
{
    /// <summary>Каталог типов фильтров</summary>
    public static IReadOnlyList<FilterKindInfo> All { get; } =
    [
        new(FilterKind.LowPass, "Фильтр нижних частот", "ФНЧ", "Пропускает частоты ниже граничной"),
        new(FilterKind.HighPass, "Фильтр верхних частот", "ФВЧ", "Пропускает частоты выше граничной"),
        new(FilterKind.BandPass, "Полосовой фильтр", "ППФ", "Пропускает заданную полосу частот"),
        new(FilterKind.BandStop, "Полосозаграждающий фильтр", "ПЗФ", "Подавляет заданную полосу частот"),
    ];

    /// <summary>Получить описание типа фильтра</summary>
    /// <param name="Kind">Тип фильтра</param>
    /// <returns>Описание типа фильтра</returns>
    public static FilterKindInfo Get(FilterKind Kind) => All.First(i => i.Kind == Kind);
}
