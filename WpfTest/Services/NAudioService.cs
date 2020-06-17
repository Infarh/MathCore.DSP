using System.Collections.Generic;
using MathCore.Annotations;
using NAudio.Wave;
using WpfTest.Models;
using WpfTest.Services.Interfaces;

namespace WpfTest.Services
{
    internal class NAudioService : IAudioService
    {
        [ItemNotNull]
        public IEnumerable<Device> GetInputDevices()
        {
            var devices_count = WaveIn.DeviceCount;
            for (var i = 0; i < devices_count; i++)
            {
                var device = WaveIn.GetCapabilities(i);
                yield return new Device(i, device.ProductName, device.Channels);
            }
        }

        [ItemNotNull]
        public IEnumerable<Device> GetOutputDevices()
        {
            var devices_count = WaveOut.DeviceCount;
            for (var i = 0; i < devices_count; i++)
            {
                var device = WaveOut.GetCapabilities(i);
                yield return new Device(i, device.ProductName, device.Channels);
            }
        }
    }
}
