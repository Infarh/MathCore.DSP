using MathCore.DI;
using MathCore.WPF.ViewModels;

using Microsoft.Extensions.DependencyInjection;

using WpfTest.ViewModels.FilterDesigners;

namespace WpfTest.ViewModels;

[Service(ServiceLifetime.Scoped)]
public class MainWindowViewModel : TitledViewModel
{
    public MainWindowViewModel() => Title = "Главное окно";

    public FilterDesign[] FilterDesigners { get; } = 
    [
        new ButterworthBandPassDesign(),
        new EllipticBandPassDesign(),
        new ChebyshevIBandPassDesign(),
        new ChebyshevIIBandPassDesign()
    ];

    #region SelectedDesigner : FilterDesign - Выбранный фильтр

    /// <summary>Выбранный фильтр</summary>
    private FilterDesign _SelectedFilterDesigner;

    /// <summary>Выбранный фильтр</summary>
    public FilterDesign SelectedFilterDesigner { get => _SelectedFilterDesigner; set => Set(ref _SelectedFilterDesigner, value); }

    #endregion
}