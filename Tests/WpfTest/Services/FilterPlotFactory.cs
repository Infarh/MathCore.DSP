#nullable enable
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Legends;
using OxyPlot.Series;

using WpfTest.Models;

namespace WpfTest.Services;

/// <summary>Цветовая палитра графиков</summary>
/// <param name="Accent">Основной цвет графиков</param>
/// <param name="Secondary">Дополнительный цвет графиков</param>
/// <param name="Warning">Цвет предупреждений</param>
/// <param name="Danger">Цвет полюсов</param>
/// <param name="Text">Цвет текста</param>
/// <param name="Muted">Цвет второстепенного текста</param>
/// <param name="Grid">Цвет линий сетки</param>
/// <param name="AxisLine">Цвет осей</param>
/// <param name="LegendBackground">Цвет фона легенды</param>
public sealed record PlotPalette(
    OxyColor Accent,
    OxyColor Secondary,
    OxyColor Warning,
    OxyColor Danger,
    OxyColor Text,
    OxyColor Muted,
    OxyColor Grid,
    OxyColor AxisLine,
    OxyColor LegendBackground);

/// <summary>Формирование и обновление графиков характеристик фильтра</summary>
public static class FilterPlotFactory
{
    /// <summary>Тёмная палитра графиков</summary>
    public static readonly PlotPalette DarkPalette = new(
        Accent: OxyColor.Parse("#5AA9FF"),
        Secondary: OxyColor.Parse("#4BE1B0"),
        Warning: OxyColor.Parse("#FFB454"),
        Danger: OxyColor.Parse("#FF6B8A"),
        Text: OxyColor.Parse("#DCE6F5"),
        Muted: OxyColor.Parse("#8894AA"),
        Grid: OxyColor.Parse("#243049"),
        AxisLine: OxyColor.Parse("#3A4A6B"),
        LegendBackground: OxyColor.Parse("#CC161C2B"));

    /// <summary>Светлая палитра графиков</summary>
    public static readonly PlotPalette LightPalette = new(
        Accent: OxyColor.Parse("#2F7DE1"),
        Secondary: OxyColor.Parse("#0E9E78"),
        Warning: OxyColor.Parse("#C77A00"),
        Danger: OxyColor.Parse("#D64560"),
        Text: OxyColor.Parse("#16202F"),
        Muted: OxyColor.Parse("#5C6B82"),
        Grid: OxyColor.Parse("#DCE3EE"),
        AxisLine: OxyColor.Parse("#B7C3D4"),
        LegendBackground: OxyColor.Parse("#E6FFFFFF"));

    /// <summary>Основной цвет графиков</summary>
    public static OxyColor Accent { get; private set; } = DarkPalette.Accent;
    /// <summary>Дополнительный цвет графиков</summary>
    public static OxyColor Secondary { get; private set; } = DarkPalette.Secondary;
    /// <summary>Цвет предупреждений</summary>
    public static OxyColor Warning { get; private set; } = DarkPalette.Warning;
    /// <summary>Цвет полюсов</summary>
    public static OxyColor Danger { get; private set; } = DarkPalette.Danger;
    /// <summary>Цвет текста</summary>
    public static OxyColor Text { get; private set; } = DarkPalette.Text;
    /// <summary>Цвет второстепенного текста</summary>
    public static OxyColor Muted { get; private set; } = DarkPalette.Muted;
    /// <summary>Цвет линий сетки</summary>
    public static OxyColor Grid { get; private set; } = DarkPalette.Grid;
    /// <summary>Цвет осей</summary>
    public static OxyColor AxisLine { get; private set; } = DarkPalette.AxisLine;

    private static PlotPalette __Palette = DarkPalette;

    private const string __MagnitudeKey = "magnitude";
    private const int __PhaseSegmentsCount = 12;
    private const string __ImpulseKey = "impulse";
    private const string __GroupDelayKey = "group-delay";
    private const string __PolesKey = "poles";
    private const string __ZerosKey = "zeros";
    private const string __WaveformKey = "waveform";

    /// <summary>Применить палитру к вновь создаваемым графикам</summary>
    /// <param name="Palette">Цветовая палитра</param>
    public static void ApplyTheme(PlotPalette Palette)
    {
        ArgumentNullException.ThrowIfNull(Palette);

        __Palette = Palette;
        Accent = Palette.Accent;
        Secondary = Palette.Secondary;
        Warning = Palette.Warning;
        Danger = Palette.Danger;
        Text = Palette.Text;
        Muted = Palette.Muted;
        Grid = Palette.Grid;
        AxisLine = Palette.AxisLine;
    }

    /// <summary>Создать график АЧХ и ФЧХ</summary>
    /// <returns>Модель графика</returns>
    public static PlotModel CreateMagnitudePhasePlot()
    {
        var model = CreateModel();

        model.Axes.Add(new LinearAxis
        {
            Position = AxisPosition.Bottom,
            Key = "freq",
            Title = "Частота, Гц",
            Minimum = 0,
            Maximum = 5_000,
            MinimumPadding = 0,
            MaximumPadding = 0,
            TextColor = Muted,
            TitleColor = Muted,
            AxislineColor = AxisLine,
            TicklineColor = AxisLine,
            MajorGridlineStyle = LineStyle.Solid,
            MajorGridlineColor = Grid,
            MinorGridlineStyle = LineStyle.None,
            FontSize = 11,
            TitleFontSize = 11,
            AxisTitleDistance = 8,
            MajorTickSize = 4,
            MinimumMinorStep = 0,
        });

        model.Axes.Add(new LinearAxis
        {
            Position = AxisPosition.Left,
            Key = "db",
            Title = "АЧХ, дБ",
            Minimum = -140,
            Maximum = 10,
            TextColor = Muted,
            TitleColor = Muted,
            AxislineColor = AxisLine,
            TicklineColor = AxisLine,
            MajorGridlineStyle = LineStyle.Solid,
            MajorGridlineColor = Grid,
            MinorGridlineStyle = LineStyle.None,
            FontSize = 11,
            TitleFontSize = 11,
            AxisTitleDistance = 8,
            MajorTickSize = 4,
        });

        model.Axes.Add(new LinearAxis
        {
            Position = AxisPosition.Right,
            Key = "deg",
            Title = "ФЧХ, °",
            Minimum = -180,
            Maximum = 180,
            MajorStep = 90,
            TextColor = Muted,
            TitleColor = Muted,
            AxislineColor = AxisLine,
            TicklineColor = AxisLine,
            MajorGridlineStyle = LineStyle.None,
            MinorGridlineStyle = LineStyle.None,
            FontSize = 11,
            TitleFontSize = 11,
            AxisTitleDistance = 8,
            MajorTickSize = 4,
        });

        model.Series.Add(new LineSeries
        {
            Title = "АЧХ",
            Tag = __MagnitudeKey,
            XAxisKey = "freq",
            YAxisKey = "db",
            Color = Accent,
            StrokeThickness = 2,
            LineJoin = LineJoin.Round,
            TrackerFormatString = "f = {2:0} Гц\n|H| = {4:0.0} дБ",
        });

        for (var i = 0; i < __PhaseSegmentsCount; i++)
            model.Series.Add(CreatePhaseSeries(i));

        AddBandEdge(model, "pass-low", "fпн");
        AddBandEdge(model, "pass-high", "fпв");
        AddBandEdge(model, "stop-low", "fзн");
        AddBandEdge(model, "stop-high", "fзв");

        return model;
    }

    /// <summary>Обновить график АЧХ и ФЧХ</summary>
    /// <param name="Model">Модель графика</param>
    /// <param name="Analysis">Результаты анализа</param>
    /// <param name="Settings">Параметры синтеза</param>
    public static void UpdateMagnitudePhasePlot(PlotModel Model, FilterAnalysis Analysis, FilterSettings Settings)
    {
        ArgumentNullException.ThrowIfNull(Model);
        ArgumentNullException.ThrowIfNull(Analysis);
        ArgumentNullException.ThrowIfNull(Settings);

        SetPoints(GetLineSeries(Model, __MagnitudeKey), Analysis.Frequencies, Analysis.MagnitudeDb);
        SetPhaseSegments(Model, Analysis.Frequencies, Analysis.PhaseDegrees);

        var frequency_axis = Model.Axes.OfType<LinearAxis>().First(a => a.Key == "freq");
        frequency_axis.Maximum = Settings.NyquistFrequency;
        frequency_axis.MajorStep = GetFrequencyStep(Settings.NyquistFrequency);

        var phase_axis = Model.Axes.OfType<LinearAxis>().First(a => a.Key == "deg");
        var finite_phase = Analysis.PhaseDegrees.Where(double.IsFinite).ToArray();
        if (finite_phase.Length > 0)
        {
            var (phase_min, phase_max) = (finite_phase.Min(), finite_phase.Max());
            var step = GetPhaseStep(phase_max - phase_min);
            phase_axis.Minimum = Math.Floor(phase_min / step) * step;
            phase_axis.Maximum = Math.Ceiling(phase_max / step) * step;
            phase_axis.MajorStep = step;
        }

        var (pass_low, pass_high, stop_low, stop_high) = GetBandEdges(Settings);
        var is_band = Settings.Kind is FilterKind.BandPass or FilterKind.BandStop;
        SetBandEdge(Model, "pass-low", pass_low ?? double.NaN, is_band ? "fпн" : "fп", Accent);
        SetBandEdge(Model, "pass-high", pass_high ?? double.NaN, "fпв", Accent);
        SetBandEdge(Model, "stop-low", stop_low ?? double.NaN, is_band ? "fзн" : "fз", Warning);
        SetBandEdge(Model, "stop-high", stop_high ?? double.NaN, "fзв", Warning);

        Model.InvalidatePlot(true);
    }

    /// <summary>Создать график импульсной характеристики</summary>
    /// <returns>Модель графика</returns>
    public static PlotModel CreateImpulsePlot()
    {
        var model = CreateModel();
        model.IsLegendVisible = false;

        model.Axes.Add(CreateValueAxis(AxisPosition.Bottom, "n", "Отсчёт", 0, 100, true));
        model.Axes.Add(CreateValueAxis(AxisPosition.Left, "amp", "Амплитуда"));

        model.Series.Add(new StemSeries
        {
            Tag = __ImpulseKey,
            XAxisKey = "n",
            YAxisKey = "amp",
            Color = Accent,
            StrokeThickness = 1.4,
            MarkerType = MarkerType.Circle,
            MarkerSize = 3,
            MarkerFill = Accent,
            MarkerStroke = Accent,
            TrackerFormatString = "n = {2:0}\nh[n] = {4:0.0000}",
        });

        return model;
    }

    /// <summary>Обновить график импульсной характеристики</summary>
    /// <param name="Model">Модель графика</param>
    /// <param name="Analysis">Результаты анализа</param>
    public static void UpdateImpulsePlot(PlotModel Model, FilterAnalysis Analysis)
    {
        ArgumentNullException.ThrowIfNull(Model);
        ArgumentNullException.ThrowIfNull(Analysis);

        var series = GetStemSeries(Model, __ImpulseKey);
        var points = series.Points;
        points.Clear();

        var impulse = Analysis.ImpulseResponse;
        for (var i = 0; i < impulse.Length; i++)
            points.Add(new DataPoint(i, impulse[i]));

        var axis = Model.Axes.OfType<LinearAxis>().First(a => a.Key == "n");
        axis.Maximum = Math.Max(16, impulse.Length - 1);

        Model.InvalidatePlot(true);
    }

    /// <summary>Создать график групповой задержки</summary>
    /// <returns>Модель графика</returns>
    public static PlotModel CreateGroupDelayPlot()
    {
        var model = CreateModel();
        model.IsLegendVisible = false;

        model.Axes.Add(new LinearAxis
        {
            Position = AxisPosition.Bottom,
            Key = "freq",
            Title = "Частота, Гц",
            Minimum = 0,
            Maximum = 5_000,
            MinimumPadding = 0,
            MaximumPadding = 0,
            TextColor = Muted,
            TitleColor = Muted,
            AxislineColor = AxisLine,
            TicklineColor = AxisLine,
            MajorGridlineStyle = LineStyle.Solid,
            MajorGridlineColor = Grid,
            FontSize = 11,
            TitleFontSize = 11,
            AxisTitleDistance = 8,
            MajorTickSize = 4,
        });

        model.Axes.Add(CreateValueAxis(AxisPosition.Left, "samples", "Задержка, отсчётов"));

        model.Series.Add(new LineSeries
        {
            Tag = __GroupDelayKey,
            XAxisKey = "freq",
            YAxisKey = "samples",
            Color = Secondary,
            StrokeThickness = 2,
            LineJoin = LineJoin.Round,
            TrackerFormatString = "f = {2:0} Гц\nτ = {4:0.00} отсч.",
        });

        return model;
    }

    /// <summary>Обновить график групповой задержки</summary>
    /// <param name="Model">Модель графика</param>
    /// <param name="Analysis">Результаты анализа</param>
    /// <param name="Settings">Параметры синтеза</param>
    public static void UpdateGroupDelayPlot(PlotModel Model, FilterAnalysis Analysis, FilterSettings Settings)
    {
        ArgumentNullException.ThrowIfNull(Model);
        ArgumentNullException.ThrowIfNull(Analysis);
        ArgumentNullException.ThrowIfNull(Settings);

        SetPoints(GetLineSeries(Model, __GroupDelayKey), Analysis.Frequencies, Analysis.GroupDelay);

        var frequency_axis = Model.Axes.OfType<LinearAxis>().First(a => a.Key == "freq");
        frequency_axis.Maximum = Settings.NyquistFrequency;
        frequency_axis.MajorStep = GetFrequencyStep(Settings.NyquistFrequency);

        Model.InvalidatePlot(true);
    }

    /// <summary>Создать график нулей и полюсов</summary>
    /// <returns>Модель графика</returns>
    public static PlotModel CreatePoleZeroPlot()
    {
        var model = CreateModel();
        model.IsLegendVisible = true;
        model.Legends[0].LegendPosition = LegendPosition.TopLeft;

        model.Axes.Add(CreateValueAxis(AxisPosition.Bottom, "re", "Re", -1.3, 1.3, true));
        model.Axes.Add(CreateValueAxis(AxisPosition.Left, "im", "Im", -1.3, 1.3));

        model.Annotations.Add(new EllipseAnnotation
        {
            X = 0,
            Y = 0,
            Width = 2,
            Height = 2,
            Stroke = AxisLine,
            StrokeThickness = 1.4,
            Fill = OxyColors.Transparent,
        });

        model.Annotations.Add(new LineAnnotation
        {
            Type = LineAnnotationType.Horizontal,
            Y = 0,
            Color = Grid,
            StrokeThickness = 1,
        });

        model.Annotations.Add(new LineAnnotation
        {
            Type = LineAnnotationType.Vertical,
            X = 0,
            Color = Grid,
            StrokeThickness = 1,
        });

        model.Series.Add(new ScatterSeries
        {
            Title = "Полюса",
            Tag = __PolesKey,
            XAxisKey = "re",
            YAxisKey = "im",
            MarkerType = MarkerType.Cross,
            MarkerSize = 9,
            MarkerStroke = Danger,
            MarkerStrokeThickness = 2.4,
            MarkerFill = OxyColors.Transparent,
            TrackerFormatString = "p = {2:0.000} {3:+0.000;-0.000}j",
        });

        model.Series.Add(new ScatterSeries
        {
            Title = "Нули",
            Tag = __ZerosKey,
            XAxisKey = "re",
            YAxisKey = "im",
            MarkerType = MarkerType.Circle,
            MarkerSize = 7,
            MarkerStroke = Secondary,
            MarkerStrokeThickness = 2.2,
            MarkerFill = OxyColors.Transparent,
            TrackerFormatString = "z = {2:0.000} {3:+0.000;-0.000}j",
        });

        return model;
    }

    /// <summary>Обновить график нулей и полюсов</summary>
    /// <param name="Model">Модель графика</param>
    /// <param name="Analysis">Результаты анализа</param>
    public static void UpdatePoleZeroPlot(PlotModel Model, FilterAnalysis Analysis)
    {
        ArgumentNullException.ThrowIfNull(Model);
        ArgumentNullException.ThrowIfNull(Analysis);

        var poles = GetScatterSeries(Model, __PolesKey);
        poles.Points.Clear();
        foreach (var pole in Analysis.Poles)
            poles.Points.Add(new ScatterPoint(pole.Real, pole.Imaginary));

        var zeros = GetScatterSeries(Model, __ZerosKey);
        zeros.Points.Clear();
        foreach (var zero in Analysis.Zeros)
            zeros.Points.Add(new ScatterPoint(zero.Real, zero.Imaginary));

        var limit = 1.3;
        foreach (var point in Analysis.Poles.Concat(Analysis.Zeros))
            limit = Math.Max(limit, Math.Max(Math.Abs(point.Real), Math.Abs(point.Imaginary)) + 0.15);

        var re_axis = Model.Axes.OfType<LinearAxis>().First(a => a.Key == "re");
        var im_axis = Model.Axes.OfType<LinearAxis>().First(a => a.Key == "im");
        re_axis.Minimum = -limit;
        re_axis.Maximum = limit;
        im_axis.Minimum = -limit;
        im_axis.Maximum = limit;

        Model.InvalidatePlot(true);
    }

    /// <summary>Создать график временной формы аудиосигнала</summary>
    /// <param name="Color">Цвет графика</param>
    /// <returns>Модель графика</returns>
    public static PlotModel CreateWaveformPlot(OxyColor Color)
    {
        var model = CreateModel();
        model.IsLegendVisible = false;

        model.Axes.Add(CreateValueAxis(AxisPosition.Bottom, "time", "Время, с", 0, 1, true));
        model.Axes.Add(CreateValueAxis(AxisPosition.Left, "amp", "Амплитуда", -1.05, 1.05));

        model.Series.Add(new AreaSeries
        {
            Tag = __WaveformKey,
            XAxisKey = "time",
            YAxisKey = "amp",
            Color = Color,
            Fill = OxyColor.FromAColor(70, Color),
            StrokeThickness = 1,
            TrackerFormatString = "t = {2:0.0000} с\na = {4:0.000}",
        });

        return model;
    }

    /// <summary>Обновить график временной формы аудиосигнала</summary>
    /// <param name="Model">Модель графика</param>
    /// <param name="Data">Аудиоданные</param>
    /// <param name="Color">Цвет графика</param>
    public static void UpdateWaveformPlot(PlotModel Model, AudioData? Data, OxyColor Color)
    {
        ArgumentNullException.ThrowIfNull(Model);

        var series = Model.Series.OfType<AreaSeries>().First(s => Equals(s.Tag, __WaveformKey));
        series.Points.Clear();
        series.Points2.Clear();

        var time_axis = Model.Axes.OfType<LinearAxis>().First(a => a.Key == "time");
        if (Data is null || Data.FrameCount == 0)
        {
            time_axis.Maximum = 1;
            Model.InvalidatePlot(true);
            return;
        }

        var (upper, lower) = BuildEnvelope(Data, 1600);
        series.Points.AddRange(upper);
        series.Points2.AddRange(lower);
        series.Color = Color;
        series.Fill = OxyColor.FromAColor(70, Color);

        time_axis.Maximum = Math.Max(0.001, Data.Duration.TotalSeconds);
        time_axis.MajorStep = GetTimeStep(Data.Duration.TotalSeconds);

        Model.InvalidatePlot(true);
    }

    /// <summary>Построить огибающую аудиосигнала для отображения</summary>
    /// <param name="Data">Аудиоданные</param>
    /// <param name="Buckets">Число столбцов огибающей</param>
    /// <returns>Верхняя и нижняя огибающие</returns>
    public static (DataPoint[] Upper, DataPoint[] Lower) BuildEnvelope(AudioData Data, int Buckets)
    {
        ArgumentNullException.ThrowIfNull(Data);
        if (Data.FrameCount == 0 || Buckets <= 0) return ([], []);

        var channels = Math.Max(1, Data.Channels);
        var frames = Data.FrameCount;
        var bucket_size = Math.Max(1, (int)Math.Ceiling((double)frames / Buckets));
        var buckets_count = (frames + bucket_size - 1) / bucket_size;

        var upper = new DataPoint[buckets_count];
        var lower = new DataPoint[buckets_count];
        var seconds_per_frame = 1d / Data.SampleRate;

        for (var bucket = 0; bucket < buckets_count; bucket++)
        {
            var start = bucket * bucket_size;
            var end = Math.Min(frames, start + bucket_size);
            var (min, max) = (double.MaxValue, double.MinValue);

            for (var frame = start; frame < end; frame++)
            {
                var value = 0d;
                for (var channel = 0; channel < channels; channel++)
                    value += Data.Samples[frame * channels + channel];
                value /= channels;

                if (value < min) min = value;
                if (value > max) max = value;
            }

            if (min > max) (min, max) = (0, 0);

            var time = (start + (end - start) * 0.5) * seconds_per_frame;
            upper[bucket] = new(time, max);
            lower[bucket] = new(time, min);
        }

        return (upper, lower);
    }

    /// <summary>Применить палитру к оформлению существующего графика</summary>
    /// <param name="Model">Модель графика</param>
    /// <param name="Palette">Цветовая палитра</param>
    public static void Restyle(PlotModel Model, PlotPalette Palette)
    {
        ArgumentNullException.ThrowIfNull(Model);
        ArgumentNullException.ThrowIfNull(Palette);

        Model.TextColor = Palette.Text;
        Model.TitleColor = Palette.Text;

        foreach (var axis in Model.Axes)
        {
            axis.TextColor = Palette.Muted;
            axis.TitleColor = Palette.Muted;
            axis.AxislineColor = Palette.AxisLine;
            axis.TicklineColor = Palette.AxisLine;
            axis.MajorGridlineColor = Palette.Grid;
            axis.MinorGridlineColor = Palette.Grid;
        }

        foreach (var legend in Model.Legends)
        {
            legend.LegendBackground = Palette.LegendBackground;
            legend.LegendBorder = Palette.Grid;
            legend.LegendTextColor = Palette.Text;
        }

        foreach (var series in Model.Series.OfType<LineSeries>())
        {
            if (series.Tag is string tag && tag.StartsWith("phase-", StringComparison.Ordinal))
                series.Color = Palette.Secondary;
            else if (Equals(series.Tag, __GroupDelayKey))
                series.Color = Palette.Secondary;
            else
                series.Color = Palette.Accent;
        }

        foreach (var series in Model.Series.OfType<StemSeries>())
        {
            series.Color = Palette.Accent;
            series.MarkerFill = Palette.Accent;
            series.MarkerStroke = Palette.Accent;
        }

        foreach (var series in Model.Series.OfType<ScatterSeries>())
        {
            if (Equals(series.Tag, __PolesKey)) series.MarkerStroke = Palette.Danger;
            else if (Equals(series.Tag, __ZerosKey)) series.MarkerStroke = Palette.Secondary;
        }

        foreach (var annotation in Model.Annotations)
        {
            annotation.TextColor = Palette.Muted;

            switch (annotation)
            {
                case EllipseAnnotation ellipse:
                    ellipse.Stroke = Palette.AxisLine;
                    break;

                case LineAnnotation { Tag: "pass-low" or "pass-high" } pass:
                    pass.Color = Palette.Accent;
                    break;

                case LineAnnotation { Tag: "stop-low" or "stop-high" } stop:
                    stop.Color = Palette.Warning;
                    break;

                case LineAnnotation line:
                    line.Color = Palette.Grid;
                    break;
            }
        }

        Model.InvalidatePlot(true);
    }

    private static PlotModel CreateModel() => new()
    {
        Background = OxyColors.Transparent,
        PlotAreaBorderColor = OxyColors.Transparent,
        PlotAreaBorderThickness = new OxyThickness(0),
        TextColor = Text,
        TitleColor = Text,
        DefaultFont = "Segoe UI",
        DefaultFontSize = 11,
        Padding = new OxyThickness(4),
        IsLegendVisible = true,
        Legends =
        {
            new Legend
            {
                LegendPosition = LegendPosition.TopRight,
                LegendPlacement = LegendPlacement.Inside,
                LegendBackground = __Palette.LegendBackground,
                LegendBorder = Grid,
                LegendBorderThickness = 1,
                LegendTextColor = Text,
                LegendFontSize = 11,
                LegendOrientation = LegendOrientation.Horizontal,
                LegendMargin = 8,
                LegendPadding = 6,
            },
        },
    };

    private static LinearAxis CreateValueAxis(
        AxisPosition Position,
        string Key,
        string Title,
        double Minimum = double.NaN,
        double Maximum = double.NaN,
        bool FixedRange = false) => new()
        {
            Position = Position,
            Key = Key,
            Title = Title,
            Minimum = Minimum,
            Maximum = Maximum,
            MinimumPadding = FixedRange ? 0 : 0.05,
            MaximumPadding = FixedRange ? 0 : 0.08,
            TextColor = Muted,
            TitleColor = Muted,
            AxislineColor = AxisLine,
            TicklineColor = AxisLine,
            MajorGridlineStyle = LineStyle.Solid,
            MajorGridlineColor = Grid,
            MinorGridlineStyle = LineStyle.None,
            FontSize = 11,
            TitleFontSize = 11,
            AxisTitleDistance = 8,
            MajorTickSize = 4,
        };

    private static void AddBandEdge(PlotModel Model, string Tag, string Text)
    {
        Model.Annotations.Add(new LineAnnotation
        {
            Type = LineAnnotationType.Vertical,
            Tag = Tag,
            XAxisKey = "freq",
            YAxisKey = "db",
            X = -1,
            Color = Accent,
            StrokeThickness = 1,
            LineStyle = LineStyle.Dash,
            Text = Text,
            TextColor = Muted,
            FontSize = 10,
            TextVerticalAlignment = VerticalAlignment.Bottom,
            TextHorizontalAlignment = HorizontalAlignment.Left,
            TextPosition = new DataPoint(-1, 0),
        });
    }

    private static void SetBandEdge(PlotModel Model, string Tag, double Frequency, string Text, OxyColor Color)
    {
        var annotation = Model.Annotations.OfType<LineAnnotation>().First(a => Equals(a.Tag, Tag));
        var is_visible = double.IsFinite(Frequency) && Frequency > 0;
        annotation.X = is_visible ? Frequency : -1;
        annotation.Text = is_visible ? Text : string.Empty;
        annotation.TextPosition = new DataPoint(annotation.X, 0);
        annotation.Color = Color;
    }

    private static (double? PassLow, double? PassHigh, double? StopLow, double? StopHigh) GetBandEdges(FilterSettings Settings)
    {
        var (pass_low, pass_high) = (Settings.PassLowFrequency, Settings.PassHighFrequency);
        var (stop_low, stop_high) = (Settings.StopLowFrequency, Settings.StopHighFrequency);

        return Settings.Kind switch
        {
            FilterKind.LowPass => (pass_low, null, stop_low, null),
            FilterKind.HighPass => (pass_low, null, stop_low, null),
            FilterKind.BandPass => (pass_low, pass_high, stop_low, stop_high),
            FilterKind.BandStop => (pass_low, pass_high, stop_low, stop_high),
            _ => (null, null, null, null),
        };
    }

    private static double GetFrequencyStep(double Nyquist)
    {
        var raw = Nyquist / 8;
        var magnitude = Math.Pow(10, Math.Floor(Math.Log10(Math.Max(raw, 1))));
        var normalized = raw / magnitude;
        var step = normalized switch
        {
            <= 1 => 1,
            <= 2 => 2,
            <= 5 => 5,
            _ => 10,
        };
        return step * magnitude;
    }

    private static double GetPhaseStep(double Range)
    {
        var raw = Math.Max(Range, 90) / 6;
        var magnitude = Math.Pow(10, Math.Floor(Math.Log10(raw)));
        var normalized = raw / magnitude;
        var step = normalized switch
        {
            <= 1 => 1,
            <= 2 => 2,
            <= 5 => 5,
            _ => 10,
        };
        return Math.Max(45, step * magnitude);
    }

    private static double GetTimeStep(double Seconds)    {
        var raw = Seconds / 8;
        if (raw <= 0) return 0.1;
        var magnitude = Math.Pow(10, Math.Floor(Math.Log10(raw)));
        var normalized = raw / magnitude;
        var step = normalized switch
        {
            <= 1 => 1,
            <= 2 => 2,
            <= 5 => 5,
            _ => 10,
        };
        return step * magnitude;
    }

    private static LineSeries CreatePhaseSeries(int Index) => new()
    {
        Title = Index == 0 ? "ФЧХ" : null,
        Tag = $"phase-{Index}",
        XAxisKey = "freq",
        YAxisKey = "deg",
        Color = Secondary,
        StrokeThickness = 1.4,
        LineJoin = LineJoin.Round,
        RenderInLegend = Index == 0,
        TrackerFormatString = "f = {2:0} Гц\nφ = {4:0.0}°",
    };

    private static void SetPhaseSegments(PlotModel Model, double[] X, double[] Y)
    {
        var segments = Model.Series
            .OfType<LineSeries>()
            .Where(s => s.Tag is string tag && tag.StartsWith("phase-", StringComparison.Ordinal))
            .ToArray();

        var segment_index = 0;
        var points = new List<DataPoint>();

        for (var i = 0; i <= Y.Length; i++)
        {
            if (i < Y.Length && double.IsFinite(Y[i]))
            {
                points.Add(new DataPoint(X[i], Y[i]));
                continue;
            }

            if (points.Count > 0 && segment_index < segments.Length)
            {
                var segment = segments[segment_index++];
                segment.Points.Clear();
                segment.Points.AddRange(points);
                segment.IsVisible = true;
            }

            points.Clear();
        }

        for (var i = segment_index; i < segments.Length; i++)
        {
            segments[i].Points.Clear();
            segments[i].IsVisible = false;
        }
    }

    private static LineSeries GetLineSeries(PlotModel Model, string Key) =>
        Model.Series.OfType<LineSeries>().First(s => Equals(s.Tag, Key));

    private static StemSeries GetStemSeries(PlotModel Model, string Key) =>
        Model.Series.OfType<StemSeries>().First(s => Equals(s.Tag, Key));

    private static ScatterSeries GetScatterSeries(PlotModel Model, string Key) =>
        Model.Series.OfType<ScatterSeries>().First(s => Equals(s.Tag, Key));

    private static void SetPoints(LineSeries Series, double[] X, double[] Y)
    {
        var points = Series.Points;
        if (points.Count != X.Length)
        {
            points.Clear();
            for (var i = 0; i < X.Length; i++) points.Add(new DataPoint(X[i], Y[i]));
        }
        else
        {
            for (var i = 0; i < X.Length; i++) points[i] = new DataPoint(X[i], Y[i]);
        }
    }
}
