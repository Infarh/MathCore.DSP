#nullable enable
using System.Windows.Input;
using System.Windows.Threading;

using MathCore.DI;
using MathCore.WPF.Commands;
using MathCore.WPF.ViewModels;

using Microsoft.Extensions.DependencyInjection;

using OxyPlot;

using WpfTest.Models;
using WpfTest.Services;
using WpfTest.Services.Interfaces;

// Аннотации обнуляемости в MathCore.WPF размечены не полностью
#pragma warning disable CS8601

namespace WpfTest.ViewModels;

/// <summary>Модель представления лаборатории исследования БИХ-фильтров</summary>
[Service(ServiceLifetime.Scoped)]
public sealed class MainWindowViewModel : TitledViewModel
{
    private const double __DefaultSamplingFrequency = 10_000;

    private readonly FilterDesignService _Design;
    private readonly FilterAnalysisService _Analysis;
    private readonly FilterCoefficientsSerializer _Serializer;
    private readonly IFileDialogService _Dialogs;
    private readonly AudioFilterService _Audio;
    private readonly DispatcherTimer _UpdateTimer;

    private FilterDesignResult? _DesignResult;
    private FilterCoefficients? _CustomCoefficients;
    private AudioData? _AudioSource;
    private AudioData? _AudioFiltered;

    private FilterKindInfo _SelectedKind = FilterKindInfo.Get(FilterKind.BandPass);
    private FilterFamilyInfo _SelectedFamily = FilterFamilyInfo.All.First(f => f.Family == FilterFamily.Butterworth);
    private IReadOnlyList<FilterFamilyInfo> _Families = [];

    private double _SamplingFrequency = __DefaultSamplingFrequency;
    private double _PassLow = 1_000;
    private double _PassHigh = 2_000;
    private double _StopLow = 500;
    private double _StopHigh = 2_500;
    private double _PassRippleDb = 1;
    private double _StopAttenuationDb = 40;
    private double _CutoffFrequency = 1_000;

    private int _FilterOrder;
    private string _FilterSummary = "—";
    private string _CoefficientsSummary = "—";
    private string _CoefficientsText = string.Empty;
    private string _StatusMessage = "Готово";
    private bool _HasError;

    private string _AudioFileName = "Файл не выбран";
    private string _AudioInfo = "Поддерживаются форматы WAV и MP3";
    private string _AudioStatus = "Загрузите аудиофайл, чтобы применить к нему текущий фильтр";
    private bool _IsAudioBusy;

    /// <summary>Инициализация модели представления лаборатории фильтров</summary>
    /// <param name="Design">Сервис синтеза фильтров</param>
    /// <param name="Analysis">Сервис расчёта характеристик</param>
    /// <param name="Serializer">Сервис сохранения коэффициентов</param>
    /// <param name="Dialogs">Сервис диалогов выбора файлов</param>
    /// <param name="Audio">Сервис обработки аудио</param>
    public MainWindowViewModel(
        FilterDesignService Design,
        FilterAnalysisService Analysis,
        FilterCoefficientsSerializer Serializer,
        IFileDialogService Dialogs,
        AudioFilterService Audio)
    {
        _Design = Design;
        _Analysis = Analysis;
        _Serializer = Serializer;
        _Dialogs = Dialogs;
        _Audio = Audio;

        Title = "Лаборатория БИХ-фильтров";

        MagnitudePhasePlot = FilterPlotFactory.CreateMagnitudePhasePlot();
        ImpulsePlot = FilterPlotFactory.CreateImpulsePlot();
        GroupDelayPlot = FilterPlotFactory.CreateGroupDelayPlot();
        PoleZeroPlot = FilterPlotFactory.CreatePoleZeroPlot();
        AudioInputPlot = FilterPlotFactory.CreateWaveformPlot(FilterPlotFactory.Muted);
        AudioOutputPlot = FilterPlotFactory.CreateWaveformPlot(FilterPlotFactory.Accent);

        _UpdateTimer = new DispatcherTimer(DispatcherPriority.Background) { Interval = TimeSpan.FromMilliseconds(30) };
        _UpdateTimer.Tick += OnUpdateTimerTick;

        SaveCoefficientsCommand = new LambdaCommand(SaveCoefficients, () => _DesignResult is not null);
        LoadCoefficientsCommand = new LambdaCommand(LoadCoefficients);
        ResetCommand = new LambdaCommand(Reset);

        OpenAudioCommand = new LambdaCommandAsync(OpenAudioAsync, () => !_IsAudioBusy);
        ApplyAudioFilterCommand = new LambdaCommandAsync(ApplyAudioFilterAsync, () => !_IsAudioBusy && _AudioSource is not null && _DesignResult is not null);
        SaveAudioCommand = new LambdaCommandAsync(SaveAudioAsync, () => !_IsAudioBusy && _AudioFiltered is not null);
        PlayOriginalCommand = new LambdaCommand(() => Play(_AudioSource), () => _AudioSource is not null);
        PlayFilteredCommand = new LambdaCommand(() => Play(_AudioFiltered), () => _AudioFiltered is not null);
        StopAudioCommand = new LambdaCommand(_Audio.Stop, () => !_IsAudioBusy);

        UpdateFamilies();
        ApplyDefaultFrequencies();
        UpdateAnalysis();
    }

    #region Графики

    /// <summary>График АЧХ и ФЧХ</summary>
    public PlotModel MagnitudePhasePlot { get; }

    /// <summary>График импульсной характеристики</summary>
    public PlotModel ImpulsePlot { get; }

    /// <summary>График групповой задержки</summary>
    public PlotModel GroupDelayPlot { get; }

    /// <summary>График нулей и полюсов</summary>
    public PlotModel PoleZeroPlot { get; }

    /// <summary>График исходного аудиосигнала</summary>
    public PlotModel AudioInputPlot { get; }

    /// <summary>График отфильтрованного аудиосигнала</summary>
    public PlotModel AudioOutputPlot { get; }

    #endregion

    #region Параметры фильтра

    /// <summary>Доступные типы фильтров</summary>
    public IReadOnlyList<FilterKindInfo> Kinds => FilterKindInfo.All;

    /// <summary>Доступные семейства фильтров для выбранного типа</summary>
    public IReadOnlyList<FilterFamilyInfo> Families
    {
        get => _Families;
        private set => Set(ref _Families, value!);
    }

    /// <summary>Выбранный тип фильтра</summary>
    public FilterKindInfo SelectedKind
    {
        get => _SelectedKind;
        set
        {
            if (value is null || ReferenceEquals(_SelectedKind, value)) return;
            _SelectedKind = value;
            OnPropertyChanged();
            OnKindChanged();
        }
    }

    /// <summary>Выбранное семейство фильтров</summary>
    public FilterFamilyInfo SelectedFamily
    {
        get => _SelectedFamily;
        set
        {
            if (value is null || ReferenceEquals(_SelectedFamily, value)) return;
            _SelectedFamily = value;
            OnPropertyChanged();
            OnFamilyChanged();
        }
    }

    /// <summary>Частота дискретизации, Гц</summary>
    public double SamplingFrequency
    {
        get => _SamplingFrequency;
        set
        {
            if (!double.IsFinite(value) || value is < 100 or > 1_000_000) return;
            if (Math.Abs(_SamplingFrequency - value) < 1e-6) return;

            var scale = value / _SamplingFrequency;
            _SamplingFrequency = value;
            OnPropertyChanged();
            OnPropertyChanged(nameof(NyquistFrequency));
            ScaleFrequencies(scale);
            ScheduleUpdate();
        }
    }

    /// <summary>Частота Найквиста, Гц</summary>
    public double NyquistFrequency => _SamplingFrequency / 2;

    /// <summary>Нижняя граница полосы пропускания, Гц</summary>
    public double PassLowFrequency
    {
        get => _PassLow;
        set { if (Set(ref _PassLow, value)) ScheduleUpdate(); }
    }

    /// <summary>Верхняя граница полосы пропускания, Гц</summary>
    public double PassHighFrequency
    {
        get => _PassHigh;
        set { if (Set(ref _PassHigh, value)) ScheduleUpdate(); }
    }

    /// <summary>Нижняя граница полосы заграждения, Гц</summary>
    public double StopLowFrequency
    {
        get => _StopLow;
        set { if (Set(ref _StopLow, value)) ScheduleUpdate(); }
    }

    /// <summary>Верхняя граница полосы заграждения, Гц</summary>
    public double StopHighFrequency
    {
        get => _StopHigh;
        set { if (Set(ref _StopHigh, value)) ScheduleUpdate(); }
    }

    /// <summary>Неравномерность АЧХ в полосе пропускания, дБ</summary>
    public double PassRippleDb
    {
        get => _PassRippleDb;
        set { if (Set(ref _PassRippleDb, value)) ScheduleUpdate(); }
    }

    /// <summary>Подавление в полосе заграждения, дБ</summary>
    public double StopAttenuationDb
    {
        get => _StopAttenuationDb;
        set { if (Set(ref _StopAttenuationDb, value)) ScheduleUpdate(); }
    }

    /// <summary>Частота среза RC-фильтра, Гц</summary>
    public double CutoffFrequency
    {
        get => _CutoffFrequency;
        set { if (Set(ref _CutoffFrequency, value)) ScheduleUpdate(); }
    }

    /// <summary>Порядок синтезированного фильтра</summary>
    public int FilterOrder
    {
        get => _FilterOrder;
        private set => Set(ref _FilterOrder, value);
    }

    /// <summary>Краткое описание синтезированного фильтра</summary>
    public string FilterSummary
    {
        get => _FilterSummary;
        private set => Set(ref _FilterSummary, value!);
    }

    /// <summary>Сводка по числу коэффициентов</summary>
    public string CoefficientsSummary
    {
        get => _CoefficientsSummary;
        private set => Set(ref _CoefficientsSummary, value!);
    }

    /// <summary>Текстовое представление коэффициентов</summary>
    public string CoefficientsText
    {
        get => _CoefficientsText;
        private set => Set(ref _CoefficientsText, value!);
    }

    /// <summary>Сообщение о состоянии расчёта</summary>
    public string StatusMessage
    {
        get => _StatusMessage;
        private set => Set(ref _StatusMessage, value!);
    }

    /// <summary>Признак ошибки синтеза фильтра</summary>
    public bool HasError
    {
        get => _HasError;
        private set => Set(ref _HasError, value);
    }

    /// <summary>Признак выбора фильтра из файла коэффициентов</summary>
    public bool IsCustomFamily => _SelectedFamily.Family == FilterFamily.Custom;

    /// <summary>Признак RC-фильтра</summary>
    public bool IsRcFamily => _SelectedFamily.Family is FilterFamily.RC or FilterFamily.RCExponential;

    /// <summary>Признак RLC-фильтра</summary>
    public bool IsRlcFamily => _SelectedFamily.Family == FilterFamily.RLC;

    /// <summary>Видимость параметров границ полос</summary>
    public bool IsEdgeParametersVisible => !IsCustomFamily && !IsRcFamily;

    /// <summary>Видимость требований к АЧХ</summary>
    public bool IsGainParametersVisible => IsEdgeParametersVisible;

    /// <summary>Видимость частоты среза RC-фильтра</summary>
    public bool IsCutoffVisible => IsRcFamily;

    /// <summary>Видимость верхней границы полосы пропускания</summary>
    public bool ShowPassHigh => _SelectedKind.Kind is FilterKind.BandPass or FilterKind.BandStop;

    /// <summary>Видимость нижней границы полосы заграждения</summary>
    public bool ShowStopLow => !IsRlcFamily;

    /// <summary>Видимость верхней границы полосы заграждения</summary>
    public bool ShowStopHigh => ShowPassHigh && !IsRlcFamily;

    /// <summary>Подпись нижней границы полосы пропускания</summary>
    public string PassLowLabel => _SelectedKind.Kind switch
    {
        FilterKind.LowPass => "Частота пропускания, Гц",
        FilterKind.HighPass => "Частота пропускания, Гц",
        FilterKind.BandPass => "Нижняя граница пропускания, Гц",
        _ => "Нижняя граница пропускания, Гц",
    };

    /// <summary>Подпись верхней границы полосы пропускания</summary>
    public string PassHighLabel => "Верхняя граница пропускания, Гц";

    /// <summary>Подпись нижней границы полосы заграждения</summary>
    public string StopLowLabel => _SelectedKind.Kind switch
    {
        FilterKind.LowPass => "Частота заграждения, Гц",
        FilterKind.HighPass => "Частота заграждения, Гц",
        FilterKind.BandPass => "Нижняя граница заграждения, Гц",
        _ => "Нижняя граница заграждения, Гц",
    };

    /// <summary>Подпись верхней границы полосы заграждения</summary>
    public string StopHighLabel => "Верхняя граница заграждения, Гц";

    #endregion

    #region Аудио

    /// <summary>Имя загруженного аудиофайла</summary>
    public string AudioFileName
    {
        get => _AudioFileName;
        private set => Set(ref _AudioFileName, value!);
    }

    /// <summary>Параметры загруженного аудиофайла</summary>
    public string AudioInfo
    {
        get => _AudioInfo;
        private set => Set(ref _AudioInfo, value!);
    }

    /// <summary>Состояние обработки аудио</summary>
    public string AudioStatus
    {
        get => _AudioStatus;
        private set => Set(ref _AudioStatus, value!);
    }

    /// <summary>Признак выполнения фоновой операции с аудио</summary>
    public bool IsAudioBusy
    {
        get => _IsAudioBusy;
        private set
        {
            if (Set(ref _IsAudioBusy, value)) CommandManager.InvalidateRequerySuggested();
        }
    }

    #endregion

    #region Команды

    /// <summary>Команда сохранения коэффициентов в файл</summary>
    public ICommand SaveCoefficientsCommand { get; }

    /// <summary>Команда загрузки коэффициентов из файла</summary>
    public ICommand LoadCoefficientsCommand { get; }

    /// <summary>Команда сброса параметров</summary>
    public ICommand ResetCommand { get; }

    /// <summary>Команда открытия аудиофайла</summary>
    public ICommand OpenAudioCommand { get; }

    /// <summary>Команда фильтрации аудиофайла</summary>
    public ICommand ApplyAudioFilterCommand { get; }

    /// <summary>Команда сохранения отфильтрованного аудиофайла</summary>
    public ICommand SaveAudioCommand { get; }

    /// <summary>Команда воспроизведения исходного сигнала</summary>
    public ICommand PlayOriginalCommand { get; }

    /// <summary>Команда воспроизведения отфильтрованного сигнала</summary>
    public ICommand PlayFilteredCommand { get; }

    /// <summary>Команда остановки воспроизведения</summary>
    public ICommand StopAudioCommand { get; }

    #endregion

    private void OnUpdateTimerTick(object? Sender, EventArgs E)
    {
        _UpdateTimer.Stop();
        UpdateAnalysis();
    }

    private void ScheduleUpdate()
    {
        _UpdateTimer.Stop();
        _UpdateTimer.Start();
    }

    private void OnKindChanged()
    {
        UpdateFamilies();

        if (!_SelectedFamily.Supports(_SelectedKind.Kind))
            SelectedFamily = _Families[0];

        ApplyDefaultFrequencies();
        NotifyLayoutChanged();
        ScheduleUpdate();
    }

    private void OnFamilyChanged()
    {
        NotifyLayoutChanged();
        ScheduleUpdate();
    }

    private void NotifyLayoutChanged()
    {
        OnPropertyChanged(nameof(IsCustomFamily));
        OnPropertyChanged(nameof(IsRcFamily));
        OnPropertyChanged(nameof(IsRlcFamily));
        OnPropertyChanged(nameof(IsEdgeParametersVisible));
        OnPropertyChanged(nameof(IsGainParametersVisible));
        OnPropertyChanged(nameof(IsCutoffVisible));
        OnPropertyChanged(nameof(ShowPassHigh));
        OnPropertyChanged(nameof(ShowStopLow));
        OnPropertyChanged(nameof(ShowStopHigh));
        OnPropertyChanged(nameof(PassLowLabel));
        OnPropertyChanged(nameof(PassHighLabel));
        OnPropertyChanged(nameof(StopLowLabel));
        OnPropertyChanged(nameof(StopHighLabel));
    }

    private void UpdateFamilies() =>
        Families = FilterFamilyInfo.All.Where(f => f.Supports(_SelectedKind.Kind)).ToArray();

    private void ApplyDefaultFrequencies()
    {
        var nyquist = NyquistFrequency;
        switch (_SelectedKind.Kind)
        {
            case FilterKind.LowPass:
                SetFrequency(ref _PassLow, 0.20 * nyquist);
                SetFrequency(ref _StopLow, 0.40 * nyquist);
                break;

            case FilterKind.HighPass:
                SetFrequency(ref _StopLow, 0.10 * nyquist);
                SetFrequency(ref _PassLow, 0.20 * nyquist);
                break;

            case FilterKind.BandPass:
                SetFrequency(ref _StopLow, 0.10 * nyquist);
                SetFrequency(ref _PassLow, 0.20 * nyquist);
                SetFrequency(ref _PassHigh, 0.40 * nyquist);
                SetFrequency(ref _StopHigh, 0.50 * nyquist);
                break;

            case FilterKind.BandStop:
                SetFrequency(ref _PassLow, 0.10 * nyquist);
                SetFrequency(ref _StopLow, 0.20 * nyquist);
                SetFrequency(ref _StopHigh, 0.40 * nyquist);
                SetFrequency(ref _PassHigh, 0.50 * nyquist);
                break;
        }

        SetFrequency(ref _CutoffFrequency, 0.20 * nyquist);
        NotifyFrequenciesChanged();

        static void SetFrequency(ref double Field, double Value) => Field = Math.Round(Value, 1);
    }

    private void ScaleFrequencies(double Scale)
    {
        if (!double.IsFinite(Scale) || Scale <= 0 || Math.Abs(Scale - 1) < 1e-9) return;
        _PassLow = Math.Round(_PassLow * Scale, 3);
        _PassHigh = Math.Round(_PassHigh * Scale, 3);
        _StopLow = Math.Round(_StopLow * Scale, 3);
        _StopHigh = Math.Round(_StopHigh * Scale, 3);
        _CutoffFrequency = Math.Round(_CutoffFrequency * Scale, 3);
        NotifyFrequenciesChanged();
    }

    private void NotifyFrequenciesChanged()
    {
        OnPropertyChanged(nameof(PassLowFrequency));
        OnPropertyChanged(nameof(PassHighFrequency));
        OnPropertyChanged(nameof(StopLowFrequency));
        OnPropertyChanged(nameof(StopHighFrequency));
        OnPropertyChanged(nameof(CutoffFrequency));
    }

    private FilterSettings CreateSettings() => new()
    {
        Kind = _SelectedKind.Kind,
        Family = _SelectedFamily.Family,
        SamplingFrequency = _SamplingFrequency,
        PassLowFrequency = _PassLow,
        PassHighFrequency = _PassHigh,
        StopLowFrequency = _StopLow,
        StopHighFrequency = _StopHigh,
        PassRippleDb = _PassRippleDb,
        StopAttenuationDb = _StopAttenuationDb,
        CutoffFrequency = _CutoffFrequency,
    };

    private void UpdateAnalysis()
    {
        try
        {
            var settings = CreateSettings();
            var result = _Design.Create(settings, _CustomCoefficients);
            var analysis = _Analysis.Calculate(result.Filter, settings.SamplingFrequency);

            _DesignResult = result;

            FilterOrder = result.Order;
            FilterSummary = result.Description;
            CoefficientsSummary = $"B: {result.Filter.B.Count} · A: {result.Filter.A.Count}";
            CoefficientsText = _Serializer.ToDisplayText(result.GetCoefficients());
            HasError = false;
            StatusMessage = $"Порядок фильтра: {result.Order}";

            FilterPlotFactory.UpdateMagnitudePhasePlot(MagnitudePhasePlot, analysis, settings);
            FilterPlotFactory.UpdateImpulsePlot(ImpulsePlot, analysis);
            FilterPlotFactory.UpdateGroupDelayPlot(GroupDelayPlot, analysis, settings);
            FilterPlotFactory.UpdatePoleZeroPlot(PoleZeroPlot, analysis);

            CommandManager.InvalidateRequerySuggested();
        }
        catch (Exception error)
        {
            HasError = true;
            StatusMessage = error.Message;
        }
    }

    private void Reset()
    {
        _SamplingFrequency = __DefaultSamplingFrequency;
        OnPropertyChanged(nameof(SamplingFrequency));
        OnPropertyChanged(nameof(NyquistFrequency));

        _SelectedKind = FilterKindInfo.Get(FilterKind.BandPass);
        OnPropertyChanged(nameof(SelectedKind));

        UpdateFamilies();
        _SelectedFamily = _Families.First(f => f.Family == FilterFamily.Butterworth);
        OnPropertyChanged(nameof(SelectedFamily));

        _PassRippleDb = 1;
        _StopAttenuationDb = 40;
        OnPropertyChanged(nameof(PassRippleDb));
        OnPropertyChanged(nameof(StopAttenuationDb));

        ApplyDefaultFrequencies();
        NotifyLayoutChanged();
        UpdateAnalysis();
    }

    private void SaveCoefficients()
    {
        if (_DesignResult is null) return;

        var coefficients = _DesignResult.GetCoefficients();
        var name = $"filter_{_SelectedKind.ShortTitle}_{_SelectedFamily.Family}_{coefficients.SamplingFrequency:F0}Hz";
        var path = _Dialogs.SaveFile(
            "Сохранить коэффициенты фильтра",
            "JSON (*.json)|*.json|Текст (*.txt)|*.txt|Все файлы (*.*)|*.*",
            name + ".json");

        if (path is null) return;

        try
        {
            _Serializer.Save(path, coefficients);
            StatusMessage = $"Коэффициенты сохранены: {System.IO.Path.GetFileName(path)}";
            HasError = false;
        }
        catch (Exception error)
        {
            HasError = true;
            StatusMessage = $"Ошибка сохранения: {error.Message}";
        }
    }

    private void LoadCoefficients()
    {
        var path = _Dialogs.OpenFile(
            "Загрузить коэффициенты фильтра",
            "Коэффициенты (*.json;*.txt)|*.json;*.txt|JSON (*.json)|*.json|Текст (*.txt)|*.txt|Все файлы (*.*)|*.*");

        if (path is null) return;

        try
        {
            var coefficients = _Serializer.Load(path);
            _CustomCoefficients = coefficients with { Description = $"Файл: {System.IO.Path.GetFileName(path)}" };

            _SelectedFamily = FilterFamilyInfo.All.First(f => f.Family == FilterFamily.Custom);
            OnPropertyChanged(nameof(SelectedFamily));
            NotifyLayoutChanged();

            if (coefficients.SamplingFrequency > 0)
            {
                if (Math.Abs(_SamplingFrequency - coefficients.SamplingFrequency) > 1e-6)
                    ScaleFrequencies(coefficients.SamplingFrequency / _SamplingFrequency);

                _SamplingFrequency = coefficients.SamplingFrequency;
                OnPropertyChanged(nameof(SamplingFrequency));
                OnPropertyChanged(nameof(NyquistFrequency));
            }

            UpdateAnalysis();
            StatusMessage = $"Загружены коэффициенты: {System.IO.Path.GetFileName(path)}";
            HasError = false;
        }
        catch (Exception error)
        {
            HasError = true;
            StatusMessage = $"Ошибка загрузки: {error.Message}";
        }
    }

    private async Task OpenAudioAsync()
    {
        var path = _Dialogs.OpenFile(
            "Открыть аудиофайл",
            "Аудио (*.wav;*.mp3)|*.wav;*.mp3|WAV (*.wav)|*.wav|MP3 (*.mp3)|*.mp3|Все файлы (*.*)|*.*");

        if (path is null) return;

        await LoadAudioAsync(path);
    }

    /// <summary>Загрузить аудиофайл по указанному пути</summary>
    /// <param name="Path">Путь к аудиофайлу</param>
    /// <returns>Задача загрузки</returns>
    public Task OpenAudioFileAsync(string Path) => LoadAudioAsync(Path);

    private async Task LoadAudioAsync(string Path)
    {
        IsAudioBusy = true;
        AudioStatus = "Загрузка аудиофайла…";

        try
        {
            var data = await Task.Run(() => _Audio.Load(Path));
            _AudioSource = data;
            _AudioFiltered = null;

            AudioFileName = System.IO.Path.GetFileName(Path);
            AudioInfo = $"{data.SampleRate} Гц · {data.Channels} кан. · {data.Duration:hh\\:mm\\:ss}";

            FilterPlotFactory.UpdateWaveformPlot(AudioInputPlot, data, FilterPlotFactory.Muted);
            FilterPlotFactory.UpdateWaveformPlot(AudioOutputPlot, null, FilterPlotFactory.Accent);

            if (Math.Abs(_SamplingFrequency - data.SampleRate) > 0.5)
            {
                var scale = data.SampleRate / _SamplingFrequency;
                _SamplingFrequency = data.SampleRate;
                OnPropertyChanged(nameof(SamplingFrequency));
                OnPropertyChanged(nameof(NyquistFrequency));
                ScaleFrequencies(scale);
                AudioStatus = $"Файл загружен. Частота дискретизации приведена к {data.SampleRate} Гц";
                UpdateAnalysis();
            }
            else
            {
                AudioStatus = "Файл загружен. Нажмите «Фильтровать»";
            }
        }
        catch (Exception error)
        {
            AudioStatus = $"Ошибка загрузки: {error.Message}";
        }
        finally
        {
            IsAudioBusy = false;
        }
    }

    private async Task ApplyAudioFilterAsync()
    {
        if (_AudioSource is null || _DesignResult is null) return;

        IsAudioBusy = true;
        AudioStatus = "Фильтрация сигнала…";

        try
        {
            var filter = _DesignResult.Filter;
            var source = _AudioSource;
            var result = await Task.Run(() => _Audio.Filter(source, filter));
            _AudioFiltered = result.Data;

            FilterPlotFactory.UpdateWaveformPlot(AudioOutputPlot, result.Data, FilterPlotFactory.Accent);

            AudioStatus = result.Gain < 1
                ? $"Готово. Применена нормировка ×{result.Gain:0.000} для предотвращения перегрузки"
                : "Готово. Нормировка не потребовалась";
        }
        catch (Exception error)
        {
            AudioStatus = $"Ошибка фильтрации: {error.Message}";
        }
        finally
        {
            IsAudioBusy = false;
        }
    }

    private async Task SaveAudioAsync()
    {
        if (_AudioFiltered is null) return;

        var name = System.IO.Path.GetFileNameWithoutExtension(_AudioFileName);
        var path = _Dialogs.SaveFile(
            "Сохранить отфильтрованный сигнал",
            "WAV (*.wav)|*.wav",
            $"filtered_{name}.wav");

        if (path is null) return;

        IsAudioBusy = true;
        AudioStatus = "Сохранение файла…";

        try
        {
            var data = _AudioFiltered;
            await Task.Run(() => _Audio.Save(path, data));
            AudioStatus = $"Сохранено: {System.IO.Path.GetFileName(path)}";
        }
        catch (Exception error)
        {
            AudioStatus = $"Ошибка сохранения: {error.Message}";
        }
        finally
        {
            IsAudioBusy = false;
        }
    }

    private void Play(AudioData? Data)
    {
        if (Data is null) return;
        try
        {
            _Audio.Play(Data);
            AudioStatus = "Воспроизведение…";
        }
        catch (Exception error)
        {
            AudioStatus = $"Ошибка воспроизведения: {error.Message}";
        }
    }

    /// <inheritdoc />
    protected override void Dispose(bool Disposing)
    {
        if (Disposing)
        {
            _UpdateTimer.Stop();
            _UpdateTimer.Tick -= OnUpdateTimerTick;
            _Audio.Stop();
        }

        base.Dispose(Disposing);
    }
}
