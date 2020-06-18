using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using MathCore.Annotations;
using MathCore.DSP.Signals;
using WpfTest.Models;

namespace WpfTest.Services.Interfaces
{
    internal interface IAudioService
    {
        [ItemNotNull]
        IEnumerable<Device> GetInputDevices();

        [ItemNotNull]
        IEnumerable<Device> GetOutputDevices();

        [ItemNotNull]
        Task<DigitalSignal> GetSignalAsync([NotNull] Device device, int SamplesCount, int SamplesRate = 44100, int BitsCount = 16, IProgress<double> Progress = null, CancellationToken Cancel = default);
    }
}
