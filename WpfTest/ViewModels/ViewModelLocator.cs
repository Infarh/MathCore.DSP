using MathCore.Hosting.WPF;

using Microsoft.Extensions.DependencyInjection;

namespace WpfTest.ViewModels
{
    class ViewModelLocator : ServiceLocatorHosted
    {
        public MicrophoneRecorderViewModel MicrophoneRecorderModel => Services.GetRequiredService<MicrophoneRecorderViewModel>();
    }
}
