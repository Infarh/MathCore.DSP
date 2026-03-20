namespace MathCore.DSP.Filters.Builders;

/// <summary>Корневой fluent-построитель IIR-фильтров</summary>
/// <param name="Dt">Период дискретизации</param>
public readonly record struct IirBuilder(double Dt)
{
    /// <summary>Построитель семейства Баттерворта</summary>
    public IirFamilyBuilder Butterworth => new(Dt, IirFilterFamily.Butterworth, ChebyshevType.I);

    /// <summary>Построитель семейства Чебышева</summary>
    /// <param name="Type">Тип фильтра Чебышева</param>
    public IirFamilyBuilder Chebyshev(ChebyshevType Type = ChebyshevType.I) => new(Dt, IirFilterFamily.Chebyshev, Type);

    /// <summary>Построитель семейства эллиптических фильтров</summary>
    public IirFamilyBuilder Elliptic => new(Dt, IirFilterFamily.Elliptic, ChebyshevType.I);

    /// <summary>Построитель RC-фильтров</summary>
    public IirRcBuilder RC => new(Dt);

    /// <summary>Построитель RLC-фильтров</summary>
    public IirRlcBuilder RLC => new(Dt);
}

/// <summary>Построитель IIR-фильтров выбранного семейства</summary>
/// <param name="Dt">Период дискретизации</param>
/// <param name="Family">Семейство фильтров</param>
/// <param name="ChebyshevType">Тип фильтра Чебышева</param>
public readonly record struct IirFamilyBuilder(double Dt, IirFilterFamily Family, ChebyshevType ChebyshevType)
{
    /// <summary>Создать построитель ФНЧ</summary>
    /// <param name="PassFrequency">Частота пропускания</param>
    /// <param name="StopFrequency">Частота заграждения</param>
    public IirSpecificationBuilder LowPass(double PassFrequency, double StopFrequency) =>
        new(Dt, Family, FrequencyPassType.LowPass, ChebyshevType, PassFrequency, StopFrequency, 0, 0, 0.891250938, 0.01);

    /// <summary>Создать построитель ФВЧ</summary>
    /// <param name="StopFrequency">Частота заграждения</param>
    /// <param name="PassFrequency">Частота пропускания</param>
    public IirSpecificationBuilder HighPass(double StopFrequency, double PassFrequency) =>
        new(Dt, Family, FrequencyPassType.HighPass, ChebyshevType, PassFrequency, StopFrequency, 0, 0, 0.891250938, 0.01);

    /// <summary>Создать построитель ППФ</summary>
    /// <param name="StopLowFrequency">Нижняя частота заграждения</param>
    /// <param name="PassLowFrequency">Нижняя частота пропускания</param>
    /// <param name="PassHighFrequency">Верхняя частота пропускания</param>
    /// <param name="StopHighFrequency">Верхняя частота заграждения</param>
    public IirSpecificationBuilder BandPass(double StopLowFrequency, double PassLowFrequency, double PassHighFrequency, double StopHighFrequency) =>
        new(Dt, Family, FrequencyPassType.BandPass, ChebyshevType, PassLowFrequency, StopLowFrequency, PassHighFrequency, StopHighFrequency, 0.891250938, 0.01);

    /// <summary>Создать построитель ПЗФ</summary>
    /// <param name="PassLowFrequency">Нижняя частота пропускания</param>
    /// <param name="StopLowFrequency">Нижняя частота заграждения</param>
    /// <param name="StopHighFrequency">Верхняя частота заграждения</param>
    /// <param name="PassHighFrequency">Верхняя частота пропускания</param>
    public IirSpecificationBuilder BandStop(double PassLowFrequency, double StopLowFrequency, double StopHighFrequency, double PassHighFrequency) =>
        new(Dt, Family, FrequencyPassType.BandStop, ChebyshevType, PassLowFrequency, StopLowFrequency, PassHighFrequency, StopHighFrequency, 0.891250938, 0.01);
}

/// <summary>Построитель IIR-фильтра по полной спецификации</summary>
/// <param name="Dt">Период дискретизации</param>
/// <param name="Family">Семейство фильтров</param>
/// <param name="PassType">Тип полосы фильтра</param>
/// <param name="ChebyshevType">Тип фильтра Чебышева</param>
/// <param name="PassLowFrequency">Нижняя или единственная частота пропускания</param>
/// <param name="StopLowFrequency">Нижняя или единственная частота заграждения</param>
/// <param name="PassHighFrequency">Верхняя частота пропускания</param>
/// <param name="StopHighFrequency">Верхняя частота заграждения</param>
/// <param name="PassGain">Коэффициент передачи в полосе пропускания</param>
/// <param name="StopGain">Коэффициент передачи в полосе заграждения</param>
public readonly record struct IirSpecificationBuilder(
    double Dt,
    IirFilterFamily Family,
    FrequencyPassType PassType,
    ChebyshevType ChebyshevType,
    double PassLowFrequency,
    double StopLowFrequency,
    double PassHighFrequency,
    double StopHighFrequency,
    double PassGain,
    double StopGain)
{
    /// <summary>Задать коэффициенты передачи в полосах</summary>
    /// <param name="PassGain">Коэффициент передачи в полосе пропускания</param>
    /// <param name="StopGain">Коэффициент передачи в полосе заграждения</param>
    public IirSpecificationBuilder WithGains(double PassGain, double StopGain) => this with { PassGain = PassGain, StopGain = StopGain };

    /// <summary>Задать тип фильтра Чебышева</summary>
    /// <param name="Type">Тип фильтра Чебышева</param>
    public IirSpecificationBuilder WithChebyshevType(ChebyshevType Type) => this with { ChebyshevType = Type };

    /// <summary>Создать экземпляр фильтра</summary>
    /// <returns>Экземпляр IIR-фильтра</returns>
    public Filter Create() => (Family, PassType) switch
    {
        (IirFilterFamily.Butterworth, FrequencyPassType.LowPass) => new ButterworthLowPass(Dt, PassLowFrequency, StopLowFrequency, PassGain, StopGain),
        (IirFilterFamily.Butterworth, FrequencyPassType.HighPass) => new ButterworthHighPass(Dt, StopLowFrequency, PassLowFrequency, PassGain, StopGain),
        (IirFilterFamily.Butterworth, FrequencyPassType.BandPass) => new ButterworthBandPass(Dt, StopLowFrequency, PassLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain),
        (IirFilterFamily.Butterworth, FrequencyPassType.BandStop) => new ButterworthBandStop(Dt, PassLowFrequency, StopLowFrequency, StopHighFrequency, PassHighFrequency, PassGain, StopGain),

        (IirFilterFamily.Chebyshev, FrequencyPassType.LowPass) => new ChebyshevLowPass(Dt, PassLowFrequency, StopLowFrequency, PassGain, StopGain, ChebyshevType),
        (IirFilterFamily.Chebyshev, FrequencyPassType.HighPass) => new ChebyshevHighPass(Dt, StopLowFrequency, PassLowFrequency, PassGain, StopGain, ChebyshevType),
        (IirFilterFamily.Chebyshev, FrequencyPassType.BandPass) => new ChebyshevBandPass(Dt, StopLowFrequency, PassLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain, ChebyshevType),
        (IirFilterFamily.Chebyshev, FrequencyPassType.BandStop) => new ChebyshevBandStop(Dt, PassLowFrequency, StopLowFrequency, StopHighFrequency, PassHighFrequency, PassGain, StopGain, ChebyshevType),

        (IirFilterFamily.Elliptic, FrequencyPassType.LowPass) => new EllipticLowPass(Dt, PassLowFrequency, StopLowFrequency, PassGain, StopGain),
        (IirFilterFamily.Elliptic, FrequencyPassType.HighPass) => new EllipticHighPass(Dt, StopLowFrequency, PassLowFrequency, PassGain, StopGain),
        (IirFilterFamily.Elliptic, FrequencyPassType.BandPass) => new EllipticBandPass(Dt, StopLowFrequency, PassLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain),
        (IirFilterFamily.Elliptic, FrequencyPassType.BandStop) => new EllipticBandStop(Dt, PassLowFrequency, StopLowFrequency, StopHighFrequency, PassHighFrequency, PassGain, StopGain),

        _ => throw new InvalidOperationException($"Неподдерживаемая комбинация Family={Family} PassType={PassType}")
    };

    /// <summary>Создать несколько независимых экземпляров фильтра</summary>
    /// <param name="Count">Количество экземпляров</param>
    /// <returns>Перечисление экземпляров фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Количество экземпляров меньше единицы</exception>
    public IEnumerable<Filter> CreateMany(int Count)
    {
        if (Count < 1)
            throw new ArgumentOutOfRangeException(nameof(Count), Count, "Количество экземпляров должно быть больше нуля");

        for (var i = 0; i < Count; i++)
            yield return Create();
    }

    /// <summary>Выполнить неявное преобразование в фильтр</summary>
    /// <param name="Builder">Построитель фильтра</param>
    public static implicit operator Filter(IirSpecificationBuilder Builder) => Builder.Create();
}

/// <summary>Построитель RC-фильтров</summary>
/// <param name="Dt">Период дискретизации</param>
public readonly record struct IirRcBuilder(double Dt)
{
    /// <summary>Создать построитель RC-ФНЧ</summary>
    /// <param name="CutoffFrequency">Частота среза</param>
    public RcLowPassSpecificationBuilder LowPass(double CutoffFrequency) => new(Dt, CutoffFrequency);

    /// <summary>Создать построитель экспоненциального RC-ФНЧ</summary>
    /// <param name="CutoffFrequency">Частота среза</param>
    public RcExponentialLowPassSpecificationBuilder LowPassExponential(double CutoffFrequency) => new(Dt, CutoffFrequency);

    /// <summary>Создать построитель RC-ФВЧ</summary>
    /// <param name="CutoffFrequency">Частота среза</param>
    public RcHighPassSpecificationBuilder HighPass(double CutoffFrequency) => new(Dt, CutoffFrequency);
}

/// <summary>Построитель RLC-фильтров</summary>
/// <param name="Dt">Период дискретизации</param>
public readonly record struct IirRlcBuilder(double Dt)
{
    /// <summary>Создать построитель RLC-ППФ</summary>
    /// <param name="ResonanceFrequency">Центральная частота</param>
    /// <param name="Bandwidth">Полоса пропускания</param>
    public RlcBandPassSpecificationBuilder BandPass(double ResonanceFrequency, double Bandwidth) => new(Dt, ResonanceFrequency, Bandwidth);

    /// <summary>Создать построитель RLC-ПЗФ</summary>
    /// <param name="ResonanceFrequency">Центральная частота</param>
    /// <param name="Bandwidth">Полоса заграждения</param>
    public RlcBandStopSpecificationBuilder BandStop(double ResonanceFrequency, double Bandwidth) => new(Dt, ResonanceFrequency, Bandwidth);
}

/// <summary>Построитель RC-ФНЧ</summary>
/// <param name="Dt">Период дискретизации</param>
/// <param name="CutoffFrequency">Частота среза</param>
public readonly record struct RcLowPassSpecificationBuilder(double Dt, double CutoffFrequency)
{
    /// <summary>Создать экземпляр фильтра</summary>
    public RCLowPass Create() => new(Dt, CutoffFrequency);

    /// <summary>Создать несколько независимых экземпляров фильтра</summary>
    /// <param name="Count">Количество экземпляров</param>
    /// <returns>Перечисление экземпляров фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Количество экземпляров меньше единицы</exception>
    public IEnumerable<RCLowPass> CreateMany(int Count)
    {
        if (Count < 1)
            throw new ArgumentOutOfRangeException(nameof(Count), Count, "Количество экземпляров должно быть больше нуля");

        for (var i = 0; i < Count; i++)
            yield return Create();
    }

    /// <summary>Выполнить неявное преобразование в фильтр</summary>
    /// <param name="Builder">Построитель фильтра</param>
    public static implicit operator RCLowPass(RcLowPassSpecificationBuilder Builder) => Builder.Create();
}

/// <summary>Построитель экспоненциального RC-ФНЧ</summary>
/// <param name="Dt">Период дискретизации</param>
/// <param name="CutoffFrequency">Частота среза</param>
public readonly record struct RcExponentialLowPassSpecificationBuilder(double Dt, double CutoffFrequency)
{
    /// <summary>Создать экземпляр фильтра</summary>
    public RCExponentialLowPass Create() => new(Dt, CutoffFrequency);

    /// <summary>Создать несколько независимых экземпляров фильтра</summary>
    /// <param name="Count">Количество экземпляров</param>
    /// <returns>Перечисление экземпляров фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Количество экземпляров меньше единицы</exception>
    public IEnumerable<RCExponentialLowPass> CreateMany(int Count)
    {
        if (Count < 1)
            throw new ArgumentOutOfRangeException(nameof(Count), Count, "Количество экземпляров должно быть больше нуля");

        for (var i = 0; i < Count; i++)
            yield return Create();
    }

    /// <summary>Выполнить неявное преобразование в фильтр</summary>
    /// <param name="Builder">Построитель фильтра</param>
    public static implicit operator RCExponentialLowPass(RcExponentialLowPassSpecificationBuilder Builder) => Builder.Create();
}

/// <summary>Построитель RC-ФВЧ</summary>
/// <param name="Dt">Период дискретизации</param>
/// <param name="CutoffFrequency">Частота среза</param>
public readonly record struct RcHighPassSpecificationBuilder(double Dt, double CutoffFrequency)
{
    /// <summary>Создать экземпляр фильтра</summary>
    public RCHighPass Create() => new(CutoffFrequency, Dt);

    /// <summary>Создать несколько независимых экземпляров фильтра</summary>
    /// <param name="Count">Количество экземпляров</param>
    /// <returns>Перечисление экземпляров фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Количество экземпляров меньше единицы</exception>
    public IEnumerable<RCHighPass> CreateMany(int Count)
    {
        if (Count < 1)
            throw new ArgumentOutOfRangeException(nameof(Count), Count, "Количество экземпляров должно быть больше нуля");

        for (var i = 0; i < Count; i++)
            yield return Create();
    }

    /// <summary>Выполнить неявное преобразование в фильтр</summary>
    /// <param name="Builder">Построитель фильтра</param>
    public static implicit operator RCHighPass(RcHighPassSpecificationBuilder Builder) => Builder.Create();
}

/// <summary>Построитель RLC-ППФ</summary>
/// <param name="Dt">Период дискретизации</param>
/// <param name="ResonanceFrequency">Центральная частота</param>
/// <param name="Bandwidth">Полоса пропускания</param>
public readonly record struct RlcBandPassSpecificationBuilder(double Dt, double ResonanceFrequency, double Bandwidth)
{
    /// <summary>Создать экземпляр фильтра</summary>
    public RLCBandPass Create() => new(ResonanceFrequency, Bandwidth, Dt);

    /// <summary>Создать несколько независимых экземпляров фильтра</summary>
    /// <param name="Count">Количество экземпляров</param>
    /// <returns>Перечисление экземпляров фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Количество экземпляров меньше единицы</exception>
    public IEnumerable<RLCBandPass> CreateMany(int Count)
    {
        if (Count < 1)
            throw new ArgumentOutOfRangeException(nameof(Count), Count, "Количество экземпляров должно быть больше нуля");

        for (var i = 0; i < Count; i++)
            yield return Create();
    }

    /// <summary>Выполнить неявное преобразование в фильтр</summary>
    /// <param name="Builder">Построитель фильтра</param>
    public static implicit operator RLCBandPass(RlcBandPassSpecificationBuilder Builder) => Builder.Create();
}

/// <summary>Построитель RLC-ПЗФ</summary>
/// <param name="Dt">Период дискретизации</param>
/// <param name="ResonanceFrequency">Центральная частота</param>
/// <param name="Bandwidth">Полоса заграждения</param>
public readonly record struct RlcBandStopSpecificationBuilder(double Dt, double ResonanceFrequency, double Bandwidth)
{
    /// <summary>Создать экземпляр фильтра</summary>
    public RLCBandStop Create() => new(ResonanceFrequency, Bandwidth, Dt);

    /// <summary>Создать несколько независимых экземпляров фильтра</summary>
    /// <param name="Count">Количество экземпляров</param>
    /// <returns>Перечисление экземпляров фильтра</returns>
    /// <exception cref="ArgumentOutOfRangeException">Количество экземпляров меньше единицы</exception>
    public IEnumerable<RLCBandStop> CreateMany(int Count)
    {
        if (Count < 1)
            throw new ArgumentOutOfRangeException(nameof(Count), Count, "Количество экземпляров должно быть больше нуля");

        for (var i = 0; i < Count; i++)
            yield return Create();
    }

    /// <summary>Выполнить неявное преобразование в фильтр</summary>
    /// <param name="Builder">Построитель фильтра</param>
    public static implicit operator RLCBandStop(RlcBandStopSpecificationBuilder Builder) => Builder.Create();
}

public enum IirFilterFamily
{
    Butterworth,
    Chebyshev,
    Elliptic
}
