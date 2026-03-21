using static System.Math;

namespace MathCore.DSP.Filters.Builders;

/// <summary>Контракт спецификации IIR-фильтра</summary>
public interface IIirFilterSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    Filter CreateFilter();
}

/// <summary>Контракт спецификации с параметрами дискретизации</summary>
public interface IIirSamplingSpecification
{
    /// <summary>Период дискретизации</summary>
    double Dt { get; }

    /// <summary>Частота дискретизации</summary>
    double Fd { get; }
}

/// <summary>Контракт спецификации с параметрами коэффициентов передачи</summary>
public interface IIirGainSpecification
{
    /// <summary>Коэффициент передачи в полосе пропускания</summary>
    double PassGain { get; }

    /// <summary>Коэффициент передачи в полосе заграждения</summary>
    double StopGain { get; }

    /// <summary>Неравномерность в полосе пропускания в децибелах</summary>
    double PassbandRippleDb { get; }

    /// <summary>Затухание в полосе заграждения в децибелах</summary>
    double StopbandAttenuationDb { get; }
}

/// <summary>Контракт спецификации с семейством фильтра</summary>
public interface IIirFamilySpecification
{
    /// <summary>Семейство фильтра</summary>
    IirFilterFamily Family { get; }
}

/// <summary>Контракт спецификации с типом полосы частот</summary>
public interface IIirPassTypeSpecification
{
    /// <summary>Тип полосы частот</summary>
    FrequencyPassType PassType { get; }
}

/// <summary>Контракт спецификации ФНЧ и ФВЧ</summary>
public interface IIirLowHighFrequenciesSpecification
{
    /// <summary>Частота пропускания</summary>
    double PassFrequency { get; }

    /// <summary>Частота заграждения</summary>
    double StopFrequency { get; }
}

/// <summary>Контракт спецификации ППФ и ПЗФ</summary>
public interface IIirBandFrequenciesSpecification
{
    /// <summary>Нижняя частота пропускания</summary>
    double PassLowFrequency { get; }

    /// <summary>Нижняя частота заграждения</summary>
    double StopLowFrequency { get; }

    /// <summary>Верхняя частота пропускания</summary>
    double PassHighFrequency { get; }

    /// <summary>Верхняя частота заграждения</summary>
    double StopHighFrequency { get; }
}

/// <summary>Контракт спецификации фильтров Чебышева</summary>
public interface IChebyshevTypeSpecification
{
    /// <summary>Тип фильтра Чебышева</summary>
    ChebyshevType ChebyshevType { get; }
}

/// <summary>Базовая спецификация IIR-фильтра</summary>
public abstract record IirFilterSpecification(
    double Dt,
    FrequencyPassType PassType,
    double PassGain,
    double StopGain)
    : IIirFilterSpecification, IIirSamplingSpecification, IIirPassTypeSpecification, IIirGainSpecification, IIirFamilySpecification
{
    /// <summary>Частота дискретизации</summary>
    public double Fd => 1 / Dt;

    /// <summary>Неравномерность в полосе пропускания в децибелах</summary>
    public double PassbandRippleDb => -20 * Log10(PassGain);

    /// <summary>Затухание в полосе заграждения в децибелах</summary>
    public double StopbandAttenuationDb => -20 * Log10(StopGain);

    /// <summary>Семейство фильтра</summary>
    public abstract IirFilterFamily Family { get; }

    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public abstract Filter CreateFilter();

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public abstract IirSpecificationBuilder ToBuilder();
}

/// <summary>Базовая спецификация семейства Баттерворта</summary>
public abstract record ButterworthSpecification(
    double Dt,
    FrequencyPassType PassType,
    double PassGain,
    double StopGain)
    : IirFilterSpecification(Dt, PassType, PassGain, StopGain)
{
    /// <summary>Семейство фильтра</summary>
    public sealed override IirFilterFamily Family => IirFilterFamily.Butterworth;
}

/// <summary>Базовая спецификация семейства Чебышева</summary>
public abstract record ChebyshevSpecification(
    double Dt,
    FrequencyPassType PassType,
    double PassGain,
    double StopGain,
    ChebyshevType ChebyshevType)
    : IirFilterSpecification(Dt, PassType, PassGain, StopGain), IChebyshevTypeSpecification
{
    /// <summary>Семейство фильтра</summary>
    public sealed override IirFilterFamily Family => IirFilterFamily.Chebyshev;
}

/// <summary>Базовая спецификация семейства эллиптических фильтров</summary>
public abstract record EllipticSpecification(
    double Dt,
    FrequencyPassType PassType,
    double PassGain,
    double StopGain)
    : IirFilterSpecification(Dt, PassType, PassGain, StopGain)
{
    /// <summary>Семейство фильтра</summary>
    public sealed override IirFilterFamily Family => IirFilterFamily.Elliptic;
}

/// <summary>Спецификация ФНЧ Баттерворта</summary>
public sealed record ButterworthLowPassSpecification(double Dt, double PassFrequency, double StopFrequency, double PassGain, double StopGain)
    : ButterworthSpecification(Dt, FrequencyPassType.LowPass, PassGain, StopGain), IIirLowHighFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new ButterworthLowPass(Dt, PassFrequency, StopFrequency, PassGain, StopGain);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType.I, PassFrequency, StopFrequency, 0, 0, PassGain, StopGain);
}

/// <summary>Спецификация ФВЧ Баттерворта</summary>
public sealed record ButterworthHighPassSpecification(double Dt, double StopFrequency, double PassFrequency, double PassGain, double StopGain)
    : ButterworthSpecification(Dt, FrequencyPassType.HighPass, PassGain, StopGain), IIirLowHighFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new ButterworthHighPass(Dt, StopFrequency, PassFrequency, PassGain, StopGain);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType.I, PassFrequency, StopFrequency, 0, 0, PassGain, StopGain);
}

/// <summary>Спецификация ППФ Баттерворта</summary>
public sealed record ButterworthBandPassSpecification(double Dt, double StopLowFrequency, double PassLowFrequency, double PassHighFrequency, double StopHighFrequency, double PassGain, double StopGain)
    : ButterworthSpecification(Dt, FrequencyPassType.BandPass, PassGain, StopGain), IIirBandFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new ButterworthBandPass(Dt, StopLowFrequency, PassLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType.I, PassLowFrequency, StopLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain);
}

/// <summary>Спецификация ПЗФ Баттерворта</summary>
public sealed record ButterworthBandStopSpecification(double Dt, double PassLowFrequency, double StopLowFrequency, double StopHighFrequency, double PassHighFrequency, double PassGain, double StopGain)
    : ButterworthSpecification(Dt, FrequencyPassType.BandStop, PassGain, StopGain), IIirBandFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new ButterworthBandStop(Dt, PassLowFrequency, StopLowFrequency, StopHighFrequency, PassHighFrequency, PassGain, StopGain);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType.I, PassLowFrequency, StopLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain);
}

/// <summary>Спецификация ФНЧ Чебышева</summary>
public sealed record ChebyshevLowPassSpecification(double Dt, double PassFrequency, double StopFrequency, double PassGain, double StopGain, ChebyshevType ChebyshevType)
    : ChebyshevSpecification(Dt, FrequencyPassType.LowPass, PassGain, StopGain, ChebyshevType), IIirLowHighFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new ChebyshevLowPass(Dt, PassFrequency, StopFrequency, PassGain, StopGain, ChebyshevType);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType, PassFrequency, StopFrequency, 0, 0, PassGain, StopGain);
}

/// <summary>Спецификация ФВЧ Чебышева</summary>
public sealed record ChebyshevHighPassSpecification(double Dt, double StopFrequency, double PassFrequency, double PassGain, double StopGain, ChebyshevType ChebyshevType)
    : ChebyshevSpecification(Dt, FrequencyPassType.HighPass, PassGain, StopGain, ChebyshevType), IIirLowHighFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new ChebyshevHighPass(Dt, StopFrequency, PassFrequency, PassGain, StopGain, ChebyshevType);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType, PassFrequency, StopFrequency, 0, 0, PassGain, StopGain);
}

/// <summary>Спецификация ППФ Чебышева</summary>
public sealed record ChebyshevBandPassSpecification(double Dt, double StopLowFrequency, double PassLowFrequency, double PassHighFrequency, double StopHighFrequency, double PassGain, double StopGain, ChebyshevType ChebyshevType)
    : ChebyshevSpecification(Dt, FrequencyPassType.BandPass, PassGain, StopGain, ChebyshevType), IIirBandFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new ChebyshevBandPass(Dt, StopLowFrequency, PassLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain, ChebyshevType);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType, PassLowFrequency, StopLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain);
}

/// <summary>Спецификация ПЗФ Чебышева</summary>
public sealed record ChebyshevBandStopSpecification(double Dt, double PassLowFrequency, double StopLowFrequency, double StopHighFrequency, double PassHighFrequency, double PassGain, double StopGain, ChebyshevType ChebyshevType)
    : ChebyshevSpecification(Dt, FrequencyPassType.BandStop, PassGain, StopGain, ChebyshevType), IIirBandFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new ChebyshevBandStop(Dt, PassLowFrequency, StopLowFrequency, StopHighFrequency, PassHighFrequency, PassGain, StopGain, ChebyshevType);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType, PassLowFrequency, StopLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain);
}

/// <summary>Спецификация ФНЧ эллиптического фильтра</summary>
public sealed record EllipticLowPassSpecification(double Dt, double PassFrequency, double StopFrequency, double PassGain, double StopGain)
    : EllipticSpecification(Dt, FrequencyPassType.LowPass, PassGain, StopGain), IIirLowHighFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new EllipticLowPass(Dt, PassFrequency, StopFrequency, PassGain, StopGain);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType.I, PassFrequency, StopFrequency, 0, 0, PassGain, StopGain);
}

/// <summary>Спецификация ФВЧ эллиптического фильтра</summary>
public sealed record EllipticHighPassSpecification(double Dt, double StopFrequency, double PassFrequency, double PassGain, double StopGain)
    : EllipticSpecification(Dt, FrequencyPassType.HighPass, PassGain, StopGain), IIirLowHighFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new EllipticHighPass(Dt, StopFrequency, PassFrequency, PassGain, StopGain);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType.I, PassFrequency, StopFrequency, 0, 0, PassGain, StopGain);
}

/// <summary>Спецификация ППФ эллиптического фильтра</summary>
public sealed record EllipticBandPassSpecification(double Dt, double StopLowFrequency, double PassLowFrequency, double PassHighFrequency, double StopHighFrequency, double PassGain, double StopGain)
    : EllipticSpecification(Dt, FrequencyPassType.BandPass, PassGain, StopGain), IIirBandFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new EllipticBandPass(Dt, StopLowFrequency, PassLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType.I, PassLowFrequency, StopLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain);
}

/// <summary>Спецификация ПЗФ эллиптического фильтра</summary>
public sealed record EllipticBandStopSpecification(double Dt, double PassLowFrequency, double StopLowFrequency, double StopHighFrequency, double PassHighFrequency, double PassGain, double StopGain)
    : EllipticSpecification(Dt, FrequencyPassType.BandStop, PassGain, StopGain), IIirBandFrequenciesSpecification
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new EllipticBandStop(Dt, PassLowFrequency, StopLowFrequency, StopHighFrequency, PassHighFrequency, PassGain, StopGain);

    /// <summary>Преобразовать спецификацию в билдер полной спецификации</summary>
    /// <returns>Билдер полной спецификации</returns>
    public override IirSpecificationBuilder ToBuilder() => new(Dt, Family, PassType, ChebyshevType.I, PassLowFrequency, StopLowFrequency, PassHighFrequency, StopHighFrequency, PassGain, StopGain);
}

/// <summary>Базовая спецификация простого IIR-фильтра на RC/RLC цепях</summary>
public abstract record SimpleIirSpecification(double Dt) : IIirFilterSpecification, IIirSamplingSpecification
{
    /// <summary>Частота дискретизации</summary>
    public double Fd => 1 / Dt;

    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public abstract Filter CreateFilter();
}

/// <summary>Спецификация RC-ФНЧ</summary>
public sealed record RcLowPassSpecification(double Dt, double CutoffFrequency) : SimpleIirSpecification(Dt)
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new RCLowPass(Dt, CutoffFrequency);
}

/// <summary>Спецификация RC-ФНЧ с экспоненциальной аппроксимацией</summary>
public sealed record RcExponentialLowPassSpecification(double Dt, double CutoffFrequency) : SimpleIirSpecification(Dt)
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new RCExponentialLowPass(Dt, CutoffFrequency);
}

/// <summary>Спецификация RC-ФВЧ</summary>
public sealed record RcHighPassSpecification(double Dt, double CutoffFrequency) : SimpleIirSpecification(Dt)
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new RCHighPass(CutoffFrequency, Dt);
}

/// <summary>Спецификация RLC-ППФ</summary>
public sealed record RlcBandPassSpecification(double Dt, double ResonanceFrequency, double Bandwidth) : SimpleIirSpecification(Dt)
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new RLCBandPass(ResonanceFrequency, Bandwidth, Dt);
}

/// <summary>Спецификация RLC-ПЗФ</summary>
public sealed record RlcBandStopSpecification(double Dt, double ResonanceFrequency, double Bandwidth) : SimpleIirSpecification(Dt)
{
    /// <summary>Создать фильтр по спецификации</summary>
    /// <returns>Экземпляр фильтра</returns>
    public override Filter CreateFilter() => new RLCBandStop(ResonanceFrequency, Bandwidth, Dt);
}
