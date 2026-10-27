#nullable enable
using MathCore.DI;
using MathCore.DSP.Filters;

using Microsoft.Extensions.DependencyInjection;

using WpfTest.Models;

namespace WpfTest.Services;

/// <summary>Синтез БИХ-фильтров по параметрам пользовательского интерфейса</summary>
[Service(ServiceLifetime.Singleton)]
public sealed class FilterDesignService
{
    /// <summary>Синтезировать фильтр по заданным параметрам</summary>
    /// <param name="Settings">Параметры синтеза</param>
    /// <param name="CustomCoefficients">Коэффициенты для семейства <see cref="FilterFamily.Custom"/></param>
    /// <returns>Результат синтеза</returns>
    /// <exception cref="FilterDesignException">Параметры фильтра некорректны</exception>
    public FilterDesignResult Create(FilterSettings Settings, FilterCoefficients? CustomCoefficients = null)
    {
        ArgumentNullException.ThrowIfNull(Settings);

        if (Settings.SamplingFrequency is < 100 or > 1_000_000)
            throw new FilterDesignException("Частота дискретизации должна быть в диапазоне 100 Гц … 1 МГц");

        var dt = 1 / Settings.SamplingFrequency;
        var gp = DbToGain(Settings.PassRippleDb);
        var gs = DbToGain(Settings.StopAttenuationDb);
        var family = FilterFamilyInfo.All.First(f => f.Family == Settings.Family);

        if (!family.Supports(Settings.Kind))
            throw new FilterDesignException($"Семейство «{family.Title}» не поддерживает тип «{FilterKindInfo.Get(Settings.Kind).Title}»");

        var description = Settings.Family == FilterFamily.Custom
            ? $"{CustomCoefficients?.Description ?? "Из файла коэффициентов"} · {FilterKindInfo.Get(Settings.Kind).ShortTitle}"
            : $"{family.Title} · {FilterKindInfo.Get(Settings.Kind).ShortTitle}";

        IIR filter = Settings.Family switch
        {
            FilterFamily.Custom => CreateCustom(CustomCoefficients),
            FilterFamily.RC => CreateRC(Settings),
            FilterFamily.RCExponential => CreateRCExponential(Settings),
            FilterFamily.RLC => CreateRLC(Settings),
            _ => CreateApproximation(Settings, dt, gp, gs),
        };

        return new(filter, Settings, description);
    }

    private static IIR CreateCustom(FilterCoefficients? Coefficients)
    {
        if (Coefficients is null)
            throw new FilterDesignException("Не загружены коэффициенты фильтра. Нажмите «Загрузить из файла…»");

        var b = Coefficients.B;
        var a = Coefficients.A;

        if (a.Length < 2)
            throw new FilterDesignException("Полином знаменателя A должен содержать не менее двух коэффициентов");

        if (b.Length == 0)
            throw new FilterDesignException("Полином числителя B не должен быть пустым");

        if (b.Length > a.Length)
            throw new FilterDesignException("Число коэффициентов B не должно превышать число коэффициентов A");

        if (a[0] == 0)
            throw new FilterDesignException("Коэффициент A[0] не должен быть равен нулю");

        if (a.Any(v => !double.IsFinite(v)) || b.Any(v => !double.IsFinite(v)))
            throw new FilterDesignException("Коэффициенты фильтра должны быть конечными числами");

        return new IIR([..b], [..a]);
    }

    private static IIR CreateRC(FilterSettings Settings)
    {
        var cutoff = Settings.CutoffFrequency;
        ValidateFrequency(cutoff, Settings.NyquistFrequency, "Частота среза RC-фильтра");

        return Settings.Kind switch
        {
            FilterKind.LowPass => new RCLowPass(1 / Settings.SamplingFrequency, cutoff),
            FilterKind.HighPass => new RCHighPass(cutoff, 1 / Settings.SamplingFrequency),
            _ => throw new FilterDesignException("RC-фильтр реализуется только для ФНЧ и ФВЧ"),
        };
    }

    private static IIR CreateRCExponential(FilterSettings Settings)
    {
        var cutoff = Settings.CutoffFrequency;
        ValidateFrequency(cutoff, Settings.NyquistFrequency, "Частота среза RC-фильтра");
        return new RCExponentialLowPass(1 / Settings.SamplingFrequency, cutoff);
    }

    private static IIR CreateRLC(FilterSettings Settings)
    {
        var (center, bandwidth) = GetResonatorParameters(Settings);
        ValidateFrequency(center, Settings.NyquistFrequency, "Центральная частота RLC-контура");

        if (bandwidth <= 0)
            throw new FilterDesignException("Полоса пропускания RLC-контура должна быть больше нуля");

        if (bandwidth >= 2 * center)
            throw new FilterDesignException("Полоса пропускания RLC-контура должна быть меньше удвоенной центральной частоты");

        var dt = 1 / Settings.SamplingFrequency;
        return Settings.Kind switch
        {
            FilterKind.BandPass => new RLCBandPass(center, bandwidth, dt),
            FilterKind.BandStop => new RLCBandStop(center, bandwidth, dt),
            _ => throw new FilterDesignException("RLC-контур реализуется только для ППФ и ПЗФ"),
        };
    }

    private static IIR CreateApproximation(FilterSettings Settings, double dt, double gp, double gs)
    {
        var (pass_low, pass_high, stop_low, stop_high) = GetBandEdges(Settings);
        var chebyshev_type = Settings.Family switch
        {
            FilterFamily.ChebyshevI => ChebyshevType.I,
            FilterFamily.ChebyshevII => ChebyshevType.II,
            FilterFamily.ChebyshevIICorrected => ChebyshevType.IICorrected,
            _ => ChebyshevType.I,
        };

        return Settings.Kind switch
        {
            FilterKind.LowPass => Settings.Family switch
            {
                FilterFamily.Butterworth => new ButterworthLowPass(dt, pass_low, stop_low, gp, gs),
                FilterFamily.Elliptic => new EllipticLowPass(dt, pass_low, stop_low, gp, gs),
                _ => new ChebyshevLowPass(dt, pass_low, stop_low, gp, gs, chebyshev_type),
            },
            FilterKind.HighPass => Settings.Family switch
            {
                FilterFamily.Butterworth => new ButterworthHighPass(dt, stop_low, pass_low, gp, gs),
                FilterFamily.Elliptic => new EllipticHighPass(dt, stop_low, pass_low, gp, gs),
                _ => new ChebyshevHighPass(dt, stop_low, pass_low, gp, gs, chebyshev_type),
            },
            FilterKind.BandPass => Settings.Family switch
            {
                FilterFamily.Butterworth => new ButterworthBandPass(dt, stop_low, pass_low, pass_high, stop_high, gp, gs),
                FilterFamily.Elliptic => new EllipticBandPass(dt, stop_low, pass_low, pass_high, stop_high, gp, gs),
                _ => new ChebyshevBandPass(dt, stop_low, pass_low, pass_high, stop_high, gp, gs, chebyshev_type),
            },
            FilterKind.BandStop => Settings.Family switch
            {
                FilterFamily.Butterworth => new ButterworthBandStop(dt, pass_low, stop_low, stop_high, pass_high, gp, gs),
                FilterFamily.Elliptic => new EllipticBandStop(dt, pass_low, stop_low, stop_high, pass_high, gp, gs),
                _ => new ChebyshevBandStop(dt, pass_low, stop_low, stop_high, pass_high, gp, gs, chebyshev_type),
            },
            _ => throw new FilterDesignException($"Неизвестный тип фильтра {Settings.Kind}"),
        };
    }

    private static (double PassLow, double PassHigh, double StopLow, double StopHigh) GetBandEdges(FilterSettings Settings)
    {
        var nyquist = Settings.NyquistFrequency;
        var (pass_low, pass_high) = (Settings.PassLowFrequency, Settings.PassHighFrequency);
        var (stop_low, stop_high) = (Settings.StopLowFrequency, Settings.StopHighFrequency);

        switch (Settings.Kind)
        {
            case FilterKind.LowPass:
                ValidateFrequency(pass_low, nyquist, "Частота пропускания");
                ValidateFrequency(stop_low, nyquist, "Частота заграждения");
                if (pass_low >= stop_low)
                    throw new FilterDesignException("Частота пропускания должна быть меньше частоты заграждения");
                break;

            case FilterKind.HighPass:
                ValidateFrequency(pass_low, nyquist, "Частота пропускания");
                ValidateFrequency(stop_low, nyquist, "Частота заграждения");
                if (stop_low >= pass_low)
                    throw new FilterDesignException("Частота заграждения должна быть меньше частоты пропускания");
                break;

            case FilterKind.BandPass:
                ValidateFrequency(stop_low, nyquist, "Нижняя частота заграждения");
                ValidateFrequency(pass_low, nyquist, "Нижняя частота пропускания");
                ValidateFrequency(pass_high, nyquist, "Верхняя частота пропускания");
                ValidateFrequency(stop_high, nyquist, "Верхняя частота заграждения");
                if (!(stop_low < pass_low && pass_low < pass_high && pass_high < stop_high))
                    throw new FilterDesignException("Частоты должны удовлетворять условию fзн < fпн < fпв < fзв");
                break;

            case FilterKind.BandStop:
                ValidateFrequency(pass_low, nyquist, "Нижняя частота пропускания");
                ValidateFrequency(stop_low, nyquist, "Нижняя частота заграждения");
                ValidateFrequency(stop_high, nyquist, "Верхняя частота заграждения");
                ValidateFrequency(pass_high, nyquist, "Верхняя частота пропускания");
                if (!(pass_low < stop_low && stop_low < stop_high && stop_high < pass_high))
                    throw new FilterDesignException("Частоты должны удовлетворять условию fпн < fзн < fзв < fпв");
                break;
        }

        return (pass_low, pass_high, stop_low, stop_high);
    }

    /// <summary>Параметры резонатора RLC по границам полосы пропускания</summary>
    /// <param name="Settings">Параметры синтеза</param>
    /// <returns>Центральная частота и полоса пропускания</returns>
    public static (double Center, double Bandwidth) GetResonatorParameters(FilterSettings Settings)
    {
        ArgumentNullException.ThrowIfNull(Settings);
        var (low, high) = (Settings.PassLowFrequency, Settings.PassHighFrequency);
        return low > 0 && high > low
            ? (Math.Sqrt(low * high), high - low)
            : (Settings.CenterFrequency, Settings.Bandwidth);
    }

    private static void ValidateFrequency(double Frequency, double Nyquist, string Name)
    {
        if (!double.IsFinite(Frequency) || Frequency <= 0 || Frequency >= Nyquist)
            throw new FilterDesignException($"{Name} должна быть в диапазоне 0 … {Nyquist:F0} Гц");
    }

    /// <summary>Перевести децибелы в коэффициент передачи</summary>
    /// <param name="Db">Значение в децибелах</param>
    /// <returns>Коэффициент передачи</returns>
    public static double DbToGain(double Db) => Math.Pow(10, -Math.Abs(Db) / 20);
}

