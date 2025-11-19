// ReSharper disable InconsistentNaming
namespace MathCore.DSP.Filters.Builders;

/// <summary>Строитель фильтров верхних частот</summary>
/// <param name="dt">Период дискретизации</param>
public readonly record struct HighPassBuilder(double dt)
{
    /// <summary>Создать строитель фильтра Баттерворта верхних частот</summary>
    /// <param name="fs">Частота заграждения</param>
    /// <param name="fp">Частота пропускания</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
    public HighPassButterworthBuilder Butterworth(double fs, double fp, double Gp, double Gs) => new(dt, fs, fp, Gp, Gs);

    /// <summary>Создать строитель фильтра Чебышева верхних частот</summary>
    /// <param name="fs">Частота заграждения</param>
    /// <param name="fp">Частота пропускания</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
    public HighPassChebyshevBuilder Chebyshev(double fs, double fp, double Gp, double Gs) => new(dt, fs, fp, Gp, Gs);
}

/// <summary>Строитель фильтра Баттерворта верхних частот</summary>
/// <param name="dt">Период дискретизации</param>
/// <param name="fs">Частота заграждения</param>
/// <param name="fp">Частота пропускания</param>
/// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
/// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
public readonly record struct HighPassButterworthBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    /// <summary>Создать фильтр Баттерворта верхних частот</summary>
    public ButterworthHighPass Create() => new(dt, fp, fs, Gp, Gs);

    /// <summary>Неявное преобразование в фильтр</summary>
    public static implicit operator Filter(HighPassButterworthBuilder builder) => builder.Create();
}

/// <summary>Строитель фильтра Чебышева верхних частот</summary>
/// <param name="dt">Период дискретизации</param>
/// <param name="fs">Частота заграждения</param>
/// <param name="fp">Частота пропускания</param>
/// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
/// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
public readonly record struct HighPassChebyshevBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    /// <summary>Создать фильтр Чебышева верхних частот</summary>
    public ChebyshevHighPass Create() => new(dt, fs, fp, Gp, Gs);

    /// <summary>Неявное преобразование в фильтр</summary>
    public static implicit operator Filter(HighPassChebyshevBuilder builder) => builder.Create();
}
