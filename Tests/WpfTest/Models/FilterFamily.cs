#nullable enable
namespace WpfTest.Models;

/// <summary>Семейство (реализация) БИХ-фильтра</summary>
public enum FilterFamily
{
    /// <summary>Фильтр Баттерворта</summary>
    Butterworth,
    /// <summary>Фильтр Чебышева первого рода</summary>
    ChebyshevI,
    /// <summary>Фильтр Чебышева второго рода</summary>
    ChebyshevII,
    /// <summary>Фильтр Чебышева второго рода с коррекцией частот</summary>
    ChebyshevIICorrected,
    /// <summary>Эллиптический фильтр Золотарева-Кауэра</summary>
    Elliptic,
    /// <summary>RC-цепочка с билинейным преобразованием</summary>
    RC,
    /// <summary>RC-цепочка с экспоненциальной аппроксимацией</summary>
    RCExponential,
    /// <summary>Резонансный RLC-контур</summary>
    RLC,
    /// <summary>Фильтр, загруженный из файла коэффициентов</summary>
    Custom,
}

/// <summary>Описание семейства фильтров для пользовательского интерфейса</summary>
/// <param name="Family">Семейство фильтров</param>
/// <param name="Title">Название</param>
/// <param name="Description">Пояснение</param>
/// <param name="Kinds">Поддерживаемые типы полос</param>
public sealed record FilterFamilyInfo(FilterFamily Family, string Title, string Description, FilterKind[] Kinds)
{
    private static readonly FilterKind[] __AllKinds =
    [
        FilterKind.LowPass,
        FilterKind.HighPass,
        FilterKind.BandPass,
        FilterKind.BandStop,
    ];

    /// <summary>Каталог семейств фильтров</summary>
    public static IReadOnlyList<FilterFamilyInfo> All { get; } =
    [
        new(FilterFamily.Butterworth, "Баттерворта",
            "Максимально плоская АЧХ в полосе пропускания", __AllKinds),
        new(FilterFamily.ChebyshevI, "Чебышева I рода",
            "Равноволновая АЧХ в полосе пропускания", __AllKinds),
        new(FilterFamily.ChebyshevII, "Чебышева II рода",
            "Равноволновое подавление в полосе заграждения", __AllKinds),
        new(FilterFamily.ChebyshevIICorrected, "Чебышева II рода (корр.)",
            "Чебышева II рода с коррекцией частотного диапазона", __AllKinds),
        new(FilterFamily.Elliptic, "Эллиптический (Золотарёва-Кауэра)",
            "Минимальный порядок при заданных требованиях к АЧХ", __AllKinds),
        new(FilterFamily.RC, "RC-цепочка",
            "Аналоговый прототип первого порядка", [FilterKind.LowPass, FilterKind.HighPass]),
        new(FilterFamily.RCExponential, "RC-экспоненциальная",
            "Дискретная RC-аппроксимация первого порядка", [FilterKind.LowPass]),
        new(FilterFamily.RLC, "RLC-контур",
            "Резонансный контур второго порядка", [FilterKind.BandPass, FilterKind.BandStop]),
        new(FilterFamily.Custom, "Из файла коэффициентов",
            "Полиномы B и A, загруженные из файла", __AllKinds),
    ];

    /// <summary>Признак поддержки указанного типа полосы</summary>
    /// <param name="Kind">Тип полосы</param>
    /// <returns>Истина, если семейство поддерживает тип полосы</returns>
    public bool Supports(FilterKind Kind) => Kinds.Contains(Kind);
}
