using System.Collections.Generic;
using WpfTest.Models;

namespace WpfTest.Services.Interfaces
{
    internal interface IAudioService
    {
        IEnumerable<Device> GetInputDevices();
        IEnumerable<Device> GetOutputDevices();
    }
}
