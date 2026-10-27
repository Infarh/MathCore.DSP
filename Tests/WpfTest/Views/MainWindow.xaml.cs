#nullable enable
using System.Runtime.InteropServices;
using System.Windows;
using System.Windows.Interop;

using WpfTest.Infrastructure;
using WpfTest.ViewModels;

namespace WpfTest.Views;

/// <summary>Главное окно лаборатории БИХ-фильтров</summary>
public partial class MainWindow : Window
{
    private const int __DwmwaUseImmersiveDarkMode = 20;
    private const int __DwmwaUseImmersiveDarkModeBefore20H1 = 19;

    private bool _IsDarkTitleBar;

    /// <summary>Инициализация главного окна</summary>
    public MainWindow()
    {
        InitializeComponent();
        Loaded += OnWindowLoaded;
        Closed += OnWindowClosed;
    }

    private void OnWindowLoaded(object Sender, RoutedEventArgs E) => ThemeManager.ThemeChanged += OnThemeChanged;

    private void OnWindowClosed(object? Sender, EventArgs E) => ThemeManager.ThemeChanged -= OnThemeChanged;

    private void OnThemeChanged(object? Sender, EventArgs E) => UpdateTitleBarTheme();

    /// <inheritdoc />
    protected override void OnSourceInitialized(EventArgs E)
    {
        base.OnSourceInitialized(E);
        UpdateTitleBarTheme();
    }

    private void UpdateTitleBarTheme()
    {
        if (!ThemeManager.IsLight == _IsDarkTitleBar) return;

        var handle = new WindowInteropHelper(this).Handle;
        if (handle == IntPtr.Zero) return;

        _IsDarkTitleBar = !ThemeManager.IsLight;
        var value = _IsDarkTitleBar ? 1 : 0;
        if (DwmSetWindowAttribute(handle, __DwmwaUseImmersiveDarkMode, ref value, sizeof(int)) != 0)
            DwmSetWindowAttribute(handle, __DwmwaUseImmersiveDarkModeBefore20H1, ref value, sizeof(int));
    }

    private void OnWindowDragOver(object Sender, DragEventArgs E)
    {
        E.Effects = GetAudioFilePath(E.Data) is null ? DragDropEffects.None : DragDropEffects.Copy;
        E.Handled = true;
    }

    private async void OnWindowDrop(object Sender, DragEventArgs E)
    {
        if (DataContext is not MainWindowViewModel view_model) return;
        if (GetAudioFilePath(E.Data) is not { } path) return;

        E.Handled = true;
        await view_model.OpenAudioFileAsync(path);
    }

    private static string? GetAudioFilePath(IDataObject Data)
    {
        if (!Data.GetDataPresent(DataFormats.FileDrop)) return null;
        if (Data.GetData(DataFormats.FileDrop) is not string[] files) return null;

        return files.FirstOrDefault(file =>
            file.EndsWith(".wav", StringComparison.OrdinalIgnoreCase)
            || file.EndsWith(".mp3", StringComparison.OrdinalIgnoreCase));
    }

    [DllImport("dwmapi.dll", PreserveSig = true)]
    private static extern int DwmSetWindowAttribute(IntPtr Hwnd, int Attribute, ref int Value, int Size);
}
