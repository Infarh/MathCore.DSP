namespace MathCore.DSP.Filters.Builders;

/// <summary>Строитель фильтров нижних частот</summary>
/// <param name="dt">Период дискретизации</param>
public readonly record struct LowPassBuilder(double dt)
{
    /// <summary>Создать строитель фильтра Баттерворта нижних частот</summary>
    /// <param name="fs">Частота заграждения</param>
    /// <param name="fp">Частота пропускания</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
    public LowPassButterworthBuilder Butterworth(double fs, double fp, double Gp, double Gs) => new(dt, fs, fp, Gp, Gs);

    /// <summary>Создать строитель фильтра Чебышева нижних частот</summary>
    /// <param name="fs">Частота заграждения</param>
    /// <param name="fp">Частота пропускания</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
    public LowPassChebyshevBuilder Chebyshev(double fs, double fp, double Gp, double Gs) => new(dt, fs, fp, Gp, Gs);

    /// <summary>Создать строитель RC-фильтра нижних частот</summary>
    /// <param name="f0">Частота среза</param>
    public LowPassRCBuilder RC(double f0) => new(dt, f0);

    /// <summary>Создать строитель RC-фильтра нижних частот с экспоненциальной аппроксимацией</summary>
    /// <param name="f0">Частота среза</param>
    public LowPassRCExponentialBuilder RCExponential(double f0) => new(dt, f0);
}

/// <summary>Строитель фильтра Баттерворта нижних частот</summary>
/// <param name="dt">Период дискретизации</param>
/// <param name="fs">Частота заграждения</param>
/// <param name="fp">Частота пропускания</param>
/// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
/// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
public readonly record struct LowPassButterworthBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    /// <summary>Создать фильтр Баттерворта нижних частот</summary>
    public ButterworthLowPass Create() => new(dt, fp, fs, Gp, Gs);

    /// <summary>Неявное преобразование в фильтр</summary>
    public static implicit operator Filter(LowPassButterworthBuilder builder) => builder.Create();
}

/// <summary>Строитель фильтра Чебышева нижних частот</summary>
/// <param name="dt">Период дискретизации</param>
/// <param name="fs">Частота заграждения</param>
/// <param name="fp">Частота пропускания</param>
/// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
/// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
public readonly record struct LowPassChebyshevBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    /// <summary>Создать фильтр Чебышева нижних частот</summary>
    public ChebyshevLowPass Create() => new(dt, fp, fs, Gp, Gs);

    /// <summary>Неявное преобразование в фильтр</summary>
    public static implicit operator Filter(LowPassChebyshevBuilder builder) => builder.Create();
}

/// <summary>Строитель эллиптического фильтра нижних частот</summary>
/// <param name="dt">Период дискретизации</param>
/// <param name="fs">Частота заграждения</param>
/// <param name="fp">Частота пропускания</param>
/// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
/// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
public readonly record struct LowPassEllipticBuilder(double dt, double fs, double fp, double Gp, double Gs)
{
    /// <summary>Создать эллиптический фильтр нижних частот</summary>
    public EllipticLowPass Create() => new(dt, fp, fs, Gp, Gs);

    /// <summary>Неявное преобразование в фильтр</summary>
    public static implicit operator Filter(LowPassEllipticBuilder builder) => builder.Create();
}

/// <summary>Строитель RC-фильтра нижних частот</summary>
/// <param name="dt">Период дискретизации</param>
/// <param name="f0">Частота среза</param>
public readonly record struct LowPassRCBuilder(double dt, double f0)
{
    /// <summary>Создать RC-фильтр нижних частот</summary>
    public RCLowPass Create() => new(dt, f0);

    /// <summary>Неявное преобразование в фильтр</summary>
    public static implicit operator Filter(LowPassRCBuilder builder) => builder.Create();
}

/// <summary>Строитель RC-фильтра нижних частот с экспоненциальной аппроксимацией</summary>
/// <param name="dt">Период дискретизации</param>
/// <param name="f0">Частота среза</param>
public readonly record struct LowPassRCExponentialBuilder(double dt, double f0)
{
    /// <summary>Создать RC-фильтр нижних частот с экспоненциальной аппроксимацией</summary>
    public RCExponentialLowPass Create() => new(dt, f0);

    /// <summary>Неявное преобразование в фильтр</summary>
    public static implicit operator Filter(LowPassRCExponentialBuilder builder) => builder.Create();
}

