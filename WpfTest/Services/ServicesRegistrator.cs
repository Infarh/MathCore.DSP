using Microsoft.Extensions.DependencyInjection;
using WpfTest.Services.Interfaces;

namespace WpfTest.Services
{
    internal static class ServicesRegistrator
    {
        public static IServiceCollection RegisterServices(this IServiceCollection services) => services
           .AddTransient<IAudioService, NAudioService>()
           .AddTransient<ISignalProcessingService, DSPService>()
        ;
    }
}
