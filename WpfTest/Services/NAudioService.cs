using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

using MathCore.Annotations;
using MathCore.DSP.Signals;

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

        public async Task<DigitalSignal> GetSignalAsync(Device device, int SamplesCount, int SamplesRate, int BitsCount, IProgress<double> Progress, CancellationToken Cancel)
        {
            using var input = new WaveIn
            {
                DeviceNumber = device.Index,
                WaveFormat = new WaveFormat(SamplesRate, BitsCount, 1)
            };

            return await input.GetSignalMono(SamplesCount, Progress, Cancel).ConfigureAwait(true);
        }
    }
}
