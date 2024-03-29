﻿using MathCore.Annotations;
using MathCore.DI;
using MathCore.DSP.Signals;

using WpfTest.Models;

namespace WpfTest.Services.Interfaces;

[Service(Implementation = typeof(NAudioService))]
internal interface IAudioService
{
    [ItemNotNull]
    IEnumerable<Device> GetInputDevices();

    [ItemNotNull]
    IEnumerable<Device> GetOutputDevices();

    [ItemNotNull]
    Task<DigitalSignal> GetSignalAsync([NotNull] Device device, int SamplesCount, int SamplesRate = 44100, int BitsCount = 16, IProgress<double> Progress = null, CancellationToken Cancel = default);

    Task PlaySignalAsync(DigitalSignal Signal, Device device, IProgress<double> Progress = null, CancellationToken Cancel = default);
}