using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;

using MathCore.DSP.Signals;

using NAudio.Wave;

using WpfTest.Models;
using WpfTest.Services.Interfaces;

namespace WpfTest.Services;

internal class NAudioService : IAudioService
{
    public IEnumerable<Device> GetInputDevices()
    {
        var devices_count = WaveIn.DeviceCount;
        for (var i = 0; i < devices_count; i++)
        {
            var device = WaveIn.GetCapabilities(i);
            yield return new Device(i, device.ProductName, device.Channels);
        }
    }

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

    public Task PlaySignalAsync(DigitalSignal Signal, Device device, IProgress<double> Progress = null, CancellationToken Cancel = default)
    {
        using var output = new WaveOut
        {
            DeviceNumber = device.Index
        };

        var tcs = new TaskCompletionSource<object>();

        output.PlaybackStopped += (s, e) => tcs.SetResult(null);

        var provider = new BufferedWaveProvider(new WaveFormat(44100, 16, 1));

        var buffer = new byte[Signal.SamplesCount * 2];


        //provider.

        output.Init(provider);

        if (Cancel.CanBeCanceled)
        {
            Cancel.Register(() => output.Stop());
            Cancel.Register(() => tcs.TrySetCanceled());
        }

        output.Play();

        return tcs.Task;
    }
}