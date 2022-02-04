using MathCore.DI;

namespace WpfTest.Services.Interfaces;

[Service(Implementation = typeof(DSPService))]
public interface ISignalProcessingService
{

}