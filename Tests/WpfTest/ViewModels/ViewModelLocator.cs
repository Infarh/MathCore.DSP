using System;

using MathCore.Hosting.WPF;

using Microsoft.Extensions.DependencyInjection;

namespace WpfTest.ViewModels;

class ViewModelLocator : ServiceLocatorHosted
{
    public MicrophoneRecorderViewModel MicrophoneRecorderModel => Services.GetRequiredService<MicrophoneRecorderViewModel>();

    public MainWindowViewModel MainModel
    {
        get
        {
            try
            {
                return Services.GetRequiredService<MainWindowViewModel>();
            }
            catch (Exception e)
            {
                LastError = e;
                throw;
            }
        }
    }

    private Exception _LastError;

    public Exception LastError
    {
        get => _LastError;
        set
        {
            if (ReferenceEquals(_LastError, value)) return;
            _LastError = value;
            OnLastErrorChanged();
        }
    }

    public event EventHandler LastErrorChanged;

    protected virtual void OnLastErrorChanged() => LastErrorChanged?.Invoke(this, EventArgs.Empty);
}