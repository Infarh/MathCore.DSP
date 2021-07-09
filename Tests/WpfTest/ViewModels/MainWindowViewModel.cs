using MathCore.Hosting;
using MathCore.WPF.ViewModels;

using Microsoft.Extensions.DependencyInjection;

using WpfTest.ViewModels.FilterDesigners;

namespace WpfTest.ViewModels
{
    [Service(ServiceLifetime.Scoped)]
    public class MainWindowViewModel : TitledViewModel
    {
        public MainWindowViewModel()
        {
            Title = "Главное окно";
        }

        public FilterDesign[] FilterDesigners { get; } = 
        {
            new ButterworthBandPassDesign(),
            new EllipticBandPassDesign(),
            new ChebyshevIBandPassDesign(),
            new ChebyshevIIBandPassDesign()
        };
    }
}
