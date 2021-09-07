using System.Linq;

using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Hosting;

namespace WpfTest
{
    public partial class App
    {
        public static ServiceDescriptor[] ServicesCollection { get; private set; }
        public static bool ServicesInitialized { get; private set; }

        protected override IHostBuilder ConfigureHostBuilderFinal(IHostBuilder builder)
        {
            var services = builder.ConfigureServices(s => ServicesCollection = s.ToArray());
            ServicesInitialized = true;
            return base.ConfigureHostBuilderFinal(builder);
        }
    }
}
