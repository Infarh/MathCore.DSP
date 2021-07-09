using MathCore.Hosting;

namespace WpfTest.Services.Interfaces
{
    [Service(Implementation = typeof(DSPService))]
    public interface ISignalProcessingService
    {

    }
}