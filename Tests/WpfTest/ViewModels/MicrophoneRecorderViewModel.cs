using System.Windows.Input;

using MathCore.Annotations;
using MathCore.DI;
using MathCore.DSP.Signals;
using MathCore.WPF;
using MathCore.WPF.Commands;
using MathCore.WPF.ViewModels;

using Microsoft.Extensions.DependencyInjection;

using WpfTest.Models;
using WpfTest.Services.Interfaces;

namespace WpfTest.ViewModels;

[Service(ServiceLifetime.Scoped)]
internal class MicrophoneRecorderViewModel(IAudioService Audio, ISignalProcessingService SignalProcessing) : ViewModel
{
    private readonly ISignalProcessingService _SignalProcessing = SignalProcessing;

    #region Title : string - Заголовок окна

    /// <summary>Заголовок окна</summary>
    private string _Title = "Заголовок главного окна";

    /// <summary>Заголовок окна</summary>
    public string Title { get => _Title; set => Set(ref _Title, value); }

    #endregion

    #region InputDevices : CollectionSelector<Device> - Записывающие устройства

    /// <summary>Записывающие устройства</summary>
    private CollectionSelector<Device> _InputDevices;

    /// <summary>Записывающие устройства</summary>
    public CollectionSelector<Device> InputDevices
    {
        get => _InputDevices;
        private set => SetValue(ref _InputDevices, value)
           .ThenIf(v => v != null, v => v.SelectFirstItem = true);
    }

    #endregion

    #region OutputDevices : CollectionSelector<Device> - Воспроизводящие устройства

    /// <summary>Воспроизводящие устройства</summary>
    private CollectionSelector<Device> _OutputDevices;

    /// <summary>Воспроизводящие устройства</summary>
    public CollectionSelector<Device> OutputDevices
    {
        get => _OutputDevices;
        private set => SetValue(ref _OutputDevices, value)
           .ThenIf(v => v != null, v => v.SelectFirstItem = true);
    }

    #endregion

    #region SignalTimeLength : double - Длина записи

    /// <summary>Длина записи</summary>
    private double _SignalTimeLength = 5;

    /// <summary>Длина записи</summary>
    public double SignalTimeLength { get => _SignalTimeLength; set => Set(ref _SignalTimeLength, value, t => t > 0); }

    #endregion

    #region RecordingProgress : double - Прогресс записи

    /// <summary>Прогресс записи</summary>
    private double _RecordingProgress;

    /// <summary>Прогресс записи</summary>
    public double RecordingProgress { get => _RecordingProgress; private set => Set(ref _RecordingProgress, value); }

    #endregion

    #region RecordedSignal : DigitalSignal - Записанный цифровой сигнал

    /// <summary>Записанный цифровой сигнал</summary>
    private DigitalSignal _RecordedSignal = MathSamplesSignal.Sin(1, 1, 0, 0.1, 10);

    /// <summary>Записанный цифровой сигнал</summary>
    public DigitalSignal RecordedSignal { get => _RecordedSignal; private set => Set(ref _RecordedSignal, value); }

    #endregion

    #region TestValue : double - Тестовое значение

    /// <summary>Тестовое значение</summary>
    private double _TestValue;

    /// <summary>Тестовое значение</summary>
    public double TestValue { get => _TestValue; set => Set(ref _TestValue, value); }

    #endregion

    #region Команды

    private ICommand _LoadInputDevicesCommand;
    [NotNull]
    public ICommand LoadInputDevicesCommand => _LoadInputDevicesCommand
        ??= new LambdaCommand(() => InputDevices = Audio.GetInputDevices().ToArray());

    private ICommand _LoadOutputDevicesCommand;
    [NotNull]
    public ICommand LoadOutputDevicesCommand => _LoadOutputDevicesCommand
        ??= new LambdaCommand(() => OutputDevices = Audio.GetOutputDevices().ToArray());

    private CancellationTokenSource _RecordingCancellation;
    private ICommand _RecordSignalCommand;
    [NotNull]
    public ICommand RecordSignalCommand => _RecordSignalCommand
        ??= new LambdaCommand(OnRecordSignalCommandExecuted, CanRecordSignalCommandExecute);

    private bool CanRecordSignalCommandExecute() => _RecordingCancellation is null;

    private async void OnRecordSignalCommandExecuted()
    {
        _RecordingCancellation = new();
        var progress = new Progress<double>(p => RecordingProgress = p);
        try
        {
            var time = _SignalTimeLength;
            const int samples_rate = 44100;
            var samples_count = (int)(time * samples_rate);
            var signal = await Audio.GetSignalAsync(InputDevices.SelectedItem, samples_count, samples_rate, Progress: progress, Cancel: _RecordingCancellation.Token);
            RecordedSignal = signal;
        }
        catch (OperationCanceledException)
        {
        }
        _RecordingCancellation = null;
        RecordingProgress = 0;
        CommandManager.InvalidateRequerySuggested();
    }

    private ICommand _RecordCancelCommand;

    public ICommand RecordCancelCommand => _RecordCancelCommand
        ??= new LambdaCommand(() => _RecordingCancellation?.Cancel(), () => _RecordingCancellation != null);

    #endregion
}