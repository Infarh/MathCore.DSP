using Microsoft.Extensions.DependencyInjection;

namespace WpfTest.ViewModels
{
    internal static class ViewModelsRegistrator
    {
        public static IServiceCollection RegisterViewModels(this IServiceCollection services) => services
           .AddSingleton<MainWindowViewModel>()
        ;
    }
}
