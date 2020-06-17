using System;
using System.Windows;
using Microsoft.Extensions.Configuration;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;
using WpfTest.ViewModels;

namespace WpfTest
{
    public partial class App
    {
        private static IHost __Host;
        public static IHost Host => __Host ??= CreateHostBuilder(Environment.GetCommandLineArgs()).Build();
        
        public static IHostBuilder CreateHostBuilder(string[] Args) =>
            Microsoft.Extensions.Hosting.Host.CreateDefaultBuilder(Args)
               .ConfigureAppConfiguration((host, cfg) => cfg.AddJsonFile("appconfig.json", true, true))
               .ConfigureServices(ConfigureServices);

        private static void ConfigureServices(HostBuilderContext host, IServiceCollection services) => services.AddSingleton<MainWindowViewModel>();

        protected override async void OnStartup(StartupEventArgs e)
        {
            var host = Host;
            base.OnStartup(e);

            await host.StartAsync().ConfigureAwait(false);
        }

        protected override async void OnExit(ExitEventArgs e)
        {
            var host = Host;
            base.OnExit(e);

            await host.StopAsync().ConfigureAwait(false);
            host.Dispose();
        }
    }
}
