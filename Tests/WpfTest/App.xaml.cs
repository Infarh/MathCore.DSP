#nullable enable
using System.Windows;
using System.Windows.Threading;

using Microsoft.Extensions.DependencyInjection;

using WpfTest.Infrastructure;
using WpfTest.Services;
using WpfTest.Services.Interfaces;

namespace WpfTest;

/// <summary>Точка входа приложения</summary>
public partial class App
{
    static App() => ConfigureServices += (_, Services) => Services.AddSingleton<IFileDialogService, FileDialogService>();

    /// <inheritdoc />
    protected override void OnStartup(StartupEventArgs E)
    {
        if (E.Args.Contains("--selftest"))
        {
            var log_path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "MathCore.DSP.WpfTest.selftest.log");
            Shutdown(SelfTest.Run(log_path));
            return;
        }

        DispatcherUnhandledException += OnDispatcherUnhandledException;
        AppDomain.CurrentDomain.UnhandledException += OnDomainUnhandledException;

        try
        {
            base.OnStartup(E);
        }
        catch (Exception error)
        {
            WriteError(error);
            throw;
        }

        // Окно создаётся после запуска хоста, когда контейнер сервисов готов
        Dispatcher.BeginInvoke(new Action(ShowMainWindow), DispatcherPriority.ApplicationIdle);
    }

    private void ShowMainWindow()
    {
        try
        {
            var locator = (ViewModels.ViewModelLocator)Resources["Locator"];
            var window = new Views.MainWindow { DataContext = locator.MainModel };
            MainWindow = window;
            window.Show();
        }
        catch (Exception error)
        {
            WriteError(error);
            Shutdown(4);
        }
    }

    private static void OnDispatcherUnhandledException(object Sender, DispatcherUnhandledExceptionEventArgs E) => WriteError(E.Exception);

    private static void OnDomainUnhandledException(object Sender, UnhandledExceptionEventArgs E) => WriteError(E.ExceptionObject as Exception);

    private static void WriteError(Exception? Error)
    {
        try
        {
            var path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "MathCore.DSP.WpfTest.error.log");
            System.IO.File.WriteAllText(path, Error?.ToString() ?? "Неизвестная ошибка");
        }
        catch
        {
            // Запись протокола ошибки не должна влиять на работу приложения
        }
    }
}
