using static System.Math;

namespace MathCore.DSP.Filters.Builders;

/// <summary>Построитель этапа задания дискретизации IIR-фильтров</summary>
public readonly record struct IirSamplingBuilder
{
    /// <summary>Задать период дискретизации</summary>
    /// <param name="SamplingPeriod">Период дискретизации</param>
    /// <returns>Построитель IIR-фильтров</returns>
    /// <exception cref="ArgumentOutOfRangeException">Период дискретизации должен быть больше нуля</exception>
    public IirBuilder WithSampling(double SamplingPeriod)
    {
        if (SamplingPeriod <= 0)
            throw new ArgumentOutOfRangeException(nameof(SamplingPeriod), SamplingPeriod, "Период дискретизации должен быть больше нуля");

        return new(SamplingPeriod);
    }

    /// <summary>Задать частоту дискретизации</summary>
    /// <param name="SamplingFrequency">Частота дискретизации</param>
    /// <returns>Построитель IIR-фильтров</returns>
    /// <exception cref="ArgumentOutOfRangeException">Частота дискретизации должна быть больше нуля</exception>
    public IirBuilder WithSamplingFrequency(double SamplingFrequency) => IirBuilder.FromSamplingFrequency(SamplingFrequency);
}

/// <summary>Корневой fluent-построитель IIR-фильтров</summary>
/// <param name="Dt">Период дискретизации</param>
public readonly record struct IirBuilder(double Dt)
{
    /// <summary>Создать построитель по частоте дискретизации</summary>
    /// <param name="SamplingFrequency">Частота дискретизации</param>
    /// <returns>Построитель IIR-фильтров</returns>
    /// <exception cref="ArgumentOutOfRangeException">Частота дискретизации должна быть больше нуля</exception>
    public static IirBuilder FromSamplingFrequency(double SamplingFrequency)
    {
        if (SamplingFrequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(SamplingFrequency), SamplingFrequency, "Частота дискретизации должна быть больше нуля");

        return new(1 / SamplingFrequency);
    }

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

    /// <summary>Задать коэффициенты передачи в полосах по неравномерности и затуханию в децибелах</summary>
    /// <param name="PassbandRippleDb">Неравномерность в полосе пропускания в дБ</param>
    /// <param name="StopbandAttenuationDb">Затухание в полосе заграждения в дБ</param>
    /// <returns>Построитель с обновлёнными коэффициентами</returns>
    /// <exception cref="ArgumentOutOfRangeException">Значения в дБ должны быть положительными</exception>
    public IirSpecificationBuilder WithGainsInDb(double PassbandRippleDb, double StopbandAttenuationDb)
    {
        if (PassbandRippleDb <= 0)
            throw new ArgumentOutOfRangeException(nameof(PassbandRippleDb), PassbandRippleDb, "Неравномерность в полосе пропускания должна быть больше нуля");

        if (StopbandAttenuationDb <= 0)
            throw new ArgumentOutOfRangeException(nameof(StopbandAttenuationDb), StopbandAttenuationDb, "Затухание в полосе заграждения должна быть больше нуля");

        var pass_gain = Pow(10, -PassbandRippleDb / 20);
        var stop_gain = Pow(10, -StopbandAttenuationDb / 20);

        return this with { PassGain = pass_gain, StopGain = stop_gain };
    }

    /// <summary>Задать коэффициент передачи в полосе пропускания по неравномерности в дБ</summary>
    /// <param name="PassbandRippleDb">Неравномерность в полосе пропускания в дБ</param>
    /// <returns>Построитель с обновлённым коэффициентом</returns>
    /// <exception cref="ArgumentOutOfRangeException">Неравномерность должна быть положительной</exception>
    public IirSpecificationBuilder WithPassbandRippleDb(double PassbandRippleDb)
    {
        if (PassbandRippleDb <= 0)
            throw new ArgumentOutOfRangeException(nameof(PassbandRippleDb), PassbandRippleDb, "Неравномерность в полосе пропускания должна быть больше нуля");

        var pass_gain = Pow(10, -PassbandRippleDb / 20);
        return this with { PassGain = pass_gain };
    }

    /// <summary>Задать коэффициент передачи в полосе заграждения по затуханию в дБ</summary>
    /// <param name="StopbandAttenuationDb">Затухание в полосе заграждения в дБ</param>
    /// <returns>Построитель с обновлённым коэффициентом</returns>
    /// <exception cref="ArgumentOutOfRangeException">Затухание должно быть положительным</exception>
    public IirSpecificationBuilder WithStopbandAttenuationDb(double StopbandAttenuationDb)
    {
        if (StopbandAttenuationDb <= 0)
            throw new ArgumentOutOfRangeException(nameof(StopbandAttenuationDb), StopbandAttenuationDb, "Затухание в полосе заграждения должно быть больше нуля");

        var stop_gain = Pow(10, -StopbandAttenuationDb / 20);
        return this with { StopGain = stop_gain };
    }

    /// <summary>Получить объект спецификации фильтра</summary>
    /// <returns>Спецификация фильтра</returns>
    public IirFilterSpecification GetSpecification() => (Family, PassType) switch
    {
        (IirFilterFamily.Butterworth, FrequencyPassType.LowPass) => new ButterworthLowPassSpecification(Dt, PassLowFrequency, StopLowFrequency, PassGain, StopGain),
        (IirFilterFamily.Butterworth, FrequencyPassType.HighPass) => new ButterworthHighPassSpecification(Dt, StopLowFrequency, PassLowFrequency, PassGain, StopGain),
        (IirFilterFamily.Butterworth, FrequencyPassType.BandPass) => new ButterworthBandPassSpecification(Dt, StopLowFrequency, PassLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain),
        (IirFilterFamily.Butterworth, FrequencyPassType.BandStop) => new ButterworthBandStopSpecification(Dt, PassLowFrequency, StopLowFrequency, StopHighFrequency, PassHighFrequency, PassGain, StopGain),

        (IirFilterFamily.Chebyshev, FrequencyPassType.LowPass) => new ChebyshevLowPassSpecification(Dt, PassLowFrequency, StopLowFrequency, PassGain, StopGain, ChebyshevType),
        (IirFilterFamily.Chebyshev, FrequencyPassType.HighPass) => new ChebyshevHighPassSpecification(Dt, StopLowFrequency, PassLowFrequency, PassGain, StopGain, ChebyshevType),
        (IirFilterFamily.Chebyshev, FrequencyPassType.BandPass) => new ChebyshevBandPassSpecification(Dt, StopLowFrequency, PassLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain, ChebyshevType),
        (IirFilterFamily.Chebyshev, FrequencyPassType.BandStop) => new ChebyshevBandStopSpecification(Dt, PassLowFrequency, StopLowFrequency, StopHighFrequency, PassHighFrequency, PassGain, StopGain, ChebyshevType),

        (IirFilterFamily.Elliptic, FrequencyPassType.LowPass) => new EllipticLowPassSpecification(Dt, PassLowFrequency, StopLowFrequency, PassGain, StopGain),
        (IirFilterFamily.Elliptic, FrequencyPassType.HighPass) => new EllipticHighPassSpecification(Dt, StopLowFrequency, PassLowFrequency, PassGain, StopGain),
        (IirFilterFamily.Elliptic, FrequencyPassType.BandPass) => new EllipticBandPassSpecification(Dt, StopLowFrequency, PassLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain),
        (IirFilterFamily.Elliptic, FrequencyPassType.BandStop) => new EllipticBandStopSpecification(Dt, PassLowFrequency, StopLowFrequency, StopHighFrequency, PassHighFrequency, PassGain, StopGain),

        _ => throw new InvalidOperationException($"Неподдерживаемая комбинация Family={Family} PassType={PassType}")
    };

    /// <summary>Создать экземпляр фильтра</summary>
    /// <returns>Экземпляр IIR-фильтра</returns>
    public Filter Create() => GetSpecification().CreateFilter();

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

    /// <summary>Попытаться создать экземпляр фильтра без выброса исключения наружу</summary>
    /// <param name="Filter">Созданный экземпляр фильтра или null</param>
    /// <returns>True, если фильтр успешно создан</returns>
    public bool TryCreate(out Filter? Filter)
    {
        try
        {
            Filter = Create();
            return true;
        }
        catch (ArgumentException)
        {
            Filter = null;
            return false;
        }
        catch (InvalidOperationException)
        {
            Filter = null;
            return false;
        }
    }

    /// <summary>Выполнить неявное преобразование в фильтр</summary>
    /// <param name="Builder">Построитель фильтра</param>
    public static implicit operator Filter(IirSpecificationBuilder Builder) => Builder.Create();
}

/// <summary>Построитель IIR-фильтров выбранного семейства</summary>
/// <param name="Dt">Период дискретизации</param>
/// <param name="Family">Семейство фильтров</param>
/// <param name="ChebyshevType">Тип фильтра Чебышева</param>
public readonly record struct IirFamilyBuilder(double Dt, IirFilterFamily Family, ChebyshevType ChebyshevType)
{
    /// <summary>Создать DSL-построитель ФНЧ</summary>
    /// <returns>DSL-построитель спецификации ФНЧ</returns>
    public IirPassbandStopbandBuilder LowPass() => new(Dt, Family, FrequencyPassType.LowPass, ChebyshevType, null, null, null, null);

    /// <summary>Создать DSL-построитель ФВЧ</summary>
    /// <returns>DSL-построитель спецификации ФВЧ</returns>
    public IirPassbandStopbandBuilder HighPass() => new(Dt, Family, FrequencyPassType.HighPass, ChebyshevType, null, null, null, null);

    /// <summary>Создать DSL-построитель ППФ</summary>
    /// <returns>DSL-построитель спецификации ППФ</returns>
    public IirPassbandStopbandBuilder BandPass() => new(Dt, Family, FrequencyPassType.BandPass, ChebyshevType, null, null, null, null);

    /// <summary>Создать DSL-построитель ПЗФ</summary>
    /// <returns>DSL-построитель спецификации ПЗФ</returns>
    public IirPassbandStopbandBuilder BandStop() => new(Dt, Family, FrequencyPassType.BandStop, ChebyshevType, null, null, null, null);

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

    /// <summary>Создать построитель ППФ по центральной частоте и ширине полосы пропускания</summary>
    /// <param name="CenterFrequency">Центральная частота полосы пропускания</param>
    /// <param name="PassbandWidth">Ширина полосы пропускания</param>
    /// <param name="TransitionWidth">Ширина переходной области с каждой стороны</param>
    /// <returns>Построитель ППФ</returns>
    /// <exception cref="ArgumentOutOfRangeException">Параметры частот должны быть положительными</exception>
    public IirSpecificationBuilder BandPassByCenter(double CenterFrequency, double PassbandWidth, double TransitionWidth)
    {
        if (CenterFrequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(CenterFrequency), CenterFrequency, "Центральная частота должна быть больше нуля");

        if (PassbandWidth <= 0)
            throw new ArgumentOutOfRangeException(nameof(PassbandWidth), PassbandWidth, "Ширина полосы пропускания должна быть больше нуля");

        if (TransitionWidth <= 0)
            throw new ArgumentOutOfRangeException(nameof(TransitionWidth), TransitionWidth, "Ширина переходной области должна быть больше нуля");

        var pass_low_frequency = CenterFrequency - PassbandWidth / 2;
        if (pass_low_frequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(PassbandWidth), PassbandWidth, "Ширина полосы пропускания приводит к неположительной нижней частоте");

        var pass_high_frequency = CenterFrequency + PassbandWidth / 2;
        var stop_low_frequency = pass_low_frequency - TransitionWidth;
        if (stop_low_frequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(TransitionWidth), TransitionWidth, "Ширина переходной области приводит к неположительной нижней частоте заграждения");

        var stop_high_frequency = pass_high_frequency + TransitionWidth;

        return BandPass(stop_low_frequency, pass_low_frequency, pass_high_frequency, stop_high_frequency);
    }

    /// <summary>Создать построитель ПЗФ по центральной частоте и ширине полосы заграждения</summary>
    /// <param name="CenterFrequency">Центральная частота полосы заграждения</param>
    /// <param name="StopbandWidth">Ширина полосы заграждения</param>
    /// <param name="TransitionWidth">Ширина переходной области с каждой стороны</param>
    /// <returns>Построитель ПЗФ</returns>
    /// <exception cref="ArgumentOutOfRangeException">Параметры частот должны быть положительными</exception>
    public IirSpecificationBuilder BandStopByCenter(double CenterFrequency, double StopbandWidth, double TransitionWidth)
    {
        if (CenterFrequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(CenterFrequency), CenterFrequency, "Центральная частота должна быть больше нуля");

        if (StopbandWidth <= 0)
            throw new ArgumentOutOfRangeException(nameof(StopbandWidth), StopbandWidth, "Ширина полосы заграждения должна быть больше нуля");

        if (TransitionWidth <= 0)
            throw new ArgumentOutOfRangeException(nameof(TransitionWidth), TransitionWidth, "Ширина переходной области должна быть больше нуля");

        var stop_low_frequency = CenterFrequency - StopbandWidth / 2;
        if (stop_low_frequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(StopbandWidth), StopbandWidth, "Ширина полосы заграждения приводит к неположительной нижней частоте");

        var stop_high_frequency = CenterFrequency + StopbandWidth / 2;
        var pass_low_frequency = stop_low_frequency - TransitionWidth;
        if (pass_low_frequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(TransitionWidth), TransitionWidth, "Ширина переходной области приводит к неположительной нижней частоте пропускания");

        var pass_high_frequency = stop_high_frequency + TransitionWidth;

        return BandStop(pass_low_frequency, stop_low_frequency, stop_high_frequency, pass_high_frequency);
    }
}

/// <summary>DSL-построитель задания полос пропускания и заграждения</summary>
/// <param name="Dt">Период дискретизации</param>
/// <param name="Family">Семейство фильтров</param>
/// <param name="PassType">Тип полосы фильтра</param>
/// <param name="ChebyshevType">Тип фильтра Чебышева</param>
/// <param name="PassLowFrequency">Нижняя или единственная частота пропускания</param>
/// <param name="StopLowFrequency">Нижняя или единственная частота заграждения</param>
/// <param name="PassHighFrequency">Верхняя частота пропускания</param>
/// <param name="StopHighFrequency">Верхняя частота заграждения</param>
public readonly record struct IirPassbandStopbandBuilder(
    double Dt,
    IirFilterFamily Family,
    FrequencyPassType PassType,
    ChebyshevType ChebyshevType,
    double? PassLowFrequency,
    double? StopLowFrequency,
    double? PassHighFrequency,
    double? StopHighFrequency)
{
    /// <summary>Задать полосу пропускания для ФНЧ и ФВЧ</summary>
    /// <param name="Frequency">Частота пропускания</param>
    /// <returns>Обновлённый DSL-построитель</returns>
    /// <exception cref="InvalidOperationException">Метод применим только к ФНЧ и ФВЧ</exception>
    /// <exception cref="ArgumentOutOfRangeException">Частота должна быть больше нуля</exception>
    public IirPassbandStopbandBuilder WithPassband(double Frequency)
    {
        if (PassType is not FrequencyPassType.LowPass and not FrequencyPassType.HighPass)
            throw new InvalidOperationException($"Метод {nameof(WithPassband)}(double) применим только к ФНЧ и ФВЧ");

        if (Frequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(Frequency), Frequency, "Частота должна быть больше нуля");

        return this with { PassLowFrequency = Frequency };
    }

    /// <summary>Задать полосу заграждения для ФНЧ и ФВЧ</summary>
    /// <param name="Frequency">Частота заграждения</param>
    /// <returns>Обновлённый DSL-построитель</returns>
    /// <exception cref="InvalidOperationException">Метод применим только к ФНЧ и ФВЧ</exception>
    /// <exception cref="ArgumentOutOfRangeException">Частота должна быть больше нуля</exception>
    public IirPassbandStopbandBuilder WithStopband(double Frequency)
    {
        if (PassType is not FrequencyPassType.LowPass and not FrequencyPassType.HighPass)
            throw new InvalidOperationException($"Метод {nameof(WithStopband)}(double) применим только к ФНЧ и ФВЧ");

        if (Frequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(Frequency), Frequency, "Частота должна быть больше нуля");

        return this with { StopLowFrequency = Frequency };
    }

    /// <summary>Задать полосу пропускания для ППФ и ПЗФ</summary>
    /// <param name="LowFrequency">Нижняя частота пропускания</param>
    /// <param name="HighFrequency">Верхняя частота пропускания</param>
    /// <returns>Обновлённый DSL-построитель</returns>
    /// <exception cref="InvalidOperationException">Метод применим только к ППФ и ПЗФ</exception>
    /// <exception cref="ArgumentOutOfRangeException">Частоты должны быть положительными и упорядоченными</exception>
    public IirPassbandStopbandBuilder WithPassband(double LowFrequency, double HighFrequency)
    {
        if (PassType is not FrequencyPassType.BandPass and not FrequencyPassType.BandStop)
            throw new InvalidOperationException($"Метод {nameof(WithPassband)}(double,double) применим только к ППФ и ПЗФ");

        if (LowFrequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(LowFrequency), LowFrequency, "Нижняя частота должна быть больше нуля");

        if (HighFrequency <= LowFrequency)
            throw new ArgumentOutOfRangeException(nameof(HighFrequency), HighFrequency, "Верхняя частота должна быть больше нижней");

        return this with { PassLowFrequency = LowFrequency, PassHighFrequency = HighFrequency };
    }

    /// <summary>Задать полосу заграждения для ППФ и ПЗФ</summary>
    /// <param name="LowFrequency">Нижняя частота заграждения</param>
    /// <param name="HighFrequency">Верхняя частота заграждения</param>
    /// <returns>Обновлённый DSL-построитель</returns>
    /// <exception cref="InvalidOperationException">Метод применим только к ППФ и ПЗФ</exception>
    /// <exception cref="ArgumentOutOfRangeException">Частоты должны быть положительными и упорядоченными</exception>
    public IirPassbandStopbandBuilder WithStopband(double LowFrequency, double HighFrequency)
    {
        if (PassType is not FrequencyPassType.BandPass and not FrequencyPassType.BandStop)
            throw new InvalidOperationException($"Метод {nameof(WithStopband)}(double,double) применим только к ППФ и ПЗФ");

        if (LowFrequency <= 0)
            throw new ArgumentOutOfRangeException(nameof(LowFrequency), LowFrequency, "Нижняя частота должна быть больше нуля");

        if (HighFrequency <= LowFrequency)
            throw new ArgumentOutOfRangeException(nameof(HighFrequency), HighFrequency, "Верхняя частота должна быть больше нижней");

        return this with { StopLowFrequency = LowFrequency, StopHighFrequency = HighFrequency };
    }

    /// <summary>Преобразовать DSL-спецификацию в построитель IIR-фильтра</summary>
    /// <returns>Построитель IIR-фильтра</returns>
    /// <exception cref="InvalidOperationException">Для выбранного типа фильтра заданы не все частоты</exception>
    public IirSpecificationBuilder ToSpecification()
    {
        static double Require(double? Value, string Name)
        {
            if (!Value.HasValue)
                throw new InvalidOperationException($"Для выбранного типа фильтра не задан параметр {Name}");

            return Value.Value;
        }

        var pass_low_frequency = Require(PassLowFrequency, nameof(PassLowFrequency));
        var stop_low_frequency = Require(StopLowFrequency, nameof(StopLowFrequency));

        return PassType switch
        {
            FrequencyPassType.LowPass or FrequencyPassType.HighPass => new(
                Dt,
                Family,
                PassType,
                ChebyshevType,
                pass_low_frequency,
                stop_low_frequency,
                0,
                0,
                0.891250938,
                0.01),

            FrequencyPassType.BandPass or FrequencyPassType.BandStop => new(
                Dt,
                Family,
                PassType,
                ChebyshevType,
                pass_low_frequency,
                stop_low_frequency,
                Require(PassHighFrequency, nameof(PassHighFrequency)),
                Require(StopHighFrequency, nameof(StopHighFrequency)),
                0.891250938,
                0.01),

            _ => throw new InvalidOperationException($"Неподдерживаемый тип фильтра {PassType}")
        };
    }

    /// <summary>Получить объект спецификации фильтра</summary>
    /// <returns>Спецификация фильтра</returns>
    public IirFilterSpecification GetSpecification() => ToSpecification().GetSpecification();

    /// <summary>Создать фильтр на основе DSL-спецификации</summary>
    /// <returns>Экземпляр IIR-фильтра</returns>
    public Filter Create() => GetSpecification().CreateFilter();
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
    /// <summary>Получить объект спецификации фильтра</summary>
    /// <returns>Спецификация фильтра</returns>
    public RcLowPassSpecification GetSpecification() => new(Dt, CutoffFrequency);

    /// <summary>Создать экземпляр фильтра</summary>
    public RCLowPass Create() => (RCLowPass)GetSpecification().CreateFilter();

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
    /// <summary>Получить объект спецификации фильтра</summary>
    /// <returns>Спецификация фильтра</returns>
    public RcExponentialLowPassSpecification GetSpecification() => new(Dt, CutoffFrequency);

    /// <summary>Создать экземпляр фильтра</summary>
    public RCExponentialLowPass Create() => (RCExponentialLowPass)GetSpecification().CreateFilter();

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
    /// <summary>Получить объект спецификации фильтра</summary>
    /// <returns>Спецификация фильтра</returns>
    public RcHighPassSpecification GetSpecification() => new(Dt, CutoffFrequency);

    /// <summary>Создать экземпляр фильтра</summary>
    public RCHighPass Create() => (RCHighPass)GetSpecification().CreateFilter();

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
    /// <summary>Получить объект спецификации фильтра</summary>
    /// <returns>Спецификация фильтра</returns>
    public RlcBandPassSpecification GetSpecification() => new(Dt, ResonanceFrequency, Bandwidth);

    /// <summary>Создать экземпляр фильтра</summary>
    public RLCBandPass Create() => (RLCBandPass)GetSpecification().CreateFilter();

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
    /// <summary>Получить объект спецификации фильтра</summary>
    /// <returns>Спецификация фильтра</returns>
    public RlcBandStopSpecification GetSpecification() => new(Dt, ResonanceFrequency, Bandwidth);

    /// <summary>Создать экземпляр фильтра</summary>
    public RLCBandStop Create() => (RLCBandStop)GetSpecification().CreateFilter();

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
