#nullable enable
using System.Windows;

using OxyPlot;

using WpfTest.Services;

namespace WpfTest.Infrastructure;

/// <summary>Управление цветовой темой приложения</summary>
public static class ThemeManager
{
    private const int __PaletteIndex = 0;

    /// <summary>Текущая тема светлая?</summary>
    public static bool IsLight { get; private set; }

    /// <summary>Событие смены темы</summary>
    public static event EventHandler? ThemeChanged;

    /// <summary>Применить тему</summary>
    /// <param name="Light">Истина для светлой темы</param>
    public static void Apply(bool Light)
    {
        IsLight = Light;
        var dictionaries = Application.Current.Resources.MergedDictionaries;
        dictionaries[__PaletteIndex] = new ResourceDictionary
        {
            Source = new Uri(Light ? "Themes/Light.xaml" : "Themes/Dark.xaml", UriKind.Relative),
        };

        FilterPlotFactory.ApplyTheme(Light ? FilterPlotFactory.LightPalette : FilterPlotFactory.DarkPalette);
        ThemeChanged?.Invoke(null, EventArgs.Empty);
    }

    /// <summary>Переключить тему</summary>
    /// <returns>Новое состояние: true — светлая тема, false — тёмная</returns>
    public static bool Toggle()
    {
        Apply(!IsLight);
        return IsLight;
    }

    /// <summary>Палитра графиков для текущей темы</summary>
    public static PlotPalette CurrentPalette => IsLight ? FilterPlotFactory.LightPalette : FilterPlotFactory.DarkPalette;
}
