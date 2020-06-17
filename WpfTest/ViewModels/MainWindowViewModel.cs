using System.Linq;
using System.Windows.Input;
using MathCore.Annotations;
using MathCore.WPF;
using MathCore.WPF.Commands;
using MathCore.WPF.ViewModels;
using WpfTest.Models;
using WpfTest.Services.Interfaces;

namespace WpfTest.ViewModels
{
    class MainWindowViewModel : ViewModel
    {
        private readonly IAudioService _Audio;
        private readonly ISignalProcessingService _SignalProcessing;

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

        #region Команды

        private ICommand _LoadInputDevicesCommand;
        [NotNull] public ICommand LoadInputDevicesCommand => _LoadInputDevicesCommand 
            ??= new LambdaCommand(() => InputDevices = _Audio.GetInputDevices().ToArray());

        private ICommand _LoadOutputDevicesCommand;
        [NotNull] public ICommand LoadOutputDevicesCommand => _LoadOutputDevicesCommand 
            ??= new LambdaCommand(() => OutputDevices = _Audio.GetOutputDevices().ToArray());

        #endregion

        public MainWindowViewModel(IAudioService Audio, ISignalProcessingService SignalProcessing)
        {
            _Audio = Audio;
            _SignalProcessing = SignalProcessing;
        }
    }
}
