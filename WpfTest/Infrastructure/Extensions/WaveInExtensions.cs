using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using MathCore.Annotations;
using MathCore.DSP.Signals;

// ReSharper disable once CheckNamespace
namespace NAudio.Wave
{
    public static class WaveInExtensions
    {
        [ItemNotNull]
        public static async Task<DigitalSignal> GetSignalMono(
            [NotNull] this WaveIn input, 
            int SamplesCount, 
            [CanBeNull] IProgress<double> Progress = null, 
            CancellationToken Cancel = default)
        {
            var last_buffer_length = -1;
            var samples_rate = input.WaveFormat.SampleRate;
            if (input.BufferMilliseconds > SamplesCount * 1000 / samples_rate)
            {
                last_buffer_length = input.BufferMilliseconds;
                input.BufferMilliseconds = SamplesCount * 1000 / samples_rate;
            }

            var tcs = new TaskCompletionSource<List<short>>();
            var samples = new List<short>(SamplesCount);
            var buffer = new short[input.BufferMilliseconds * samples_rate / 1000];

            input.DataAvailable += (s, e) =>
            {
                var recorded = e.BytesRecorded / 2;
                var current_count = samples.Count;
                var capacity = samples.Capacity;
                var necessary_count = capacity - current_count;
                Progress?.Report((double)current_count / capacity);
                if (necessary_count == 0)
                {
                    ((WaveIn)s)?.StopRecording();
                    return;
                }
                if (necessary_count >= recorded && recorded == buffer.Length)
                {
                    Buffer.BlockCopy(e.Buffer, 0, buffer, 0, e.BytesRecorded);
                    samples.AddRange(buffer);
                    return;
                }

                Buffer.BlockCopy(e.Buffer, 0, buffer, 0, Math.Min(necessary_count * 2, e.Buffer.Length));
                for (var i = 0; i < necessary_count && i < buffer.Length; i++)
                    samples.Add(buffer[i]);
                ((WaveIn)s)?.StopRecording();
            };
            input.RecordingStopped += (s, e) => tcs.TrySetResult(samples);

            if (Cancel.CanBeCanceled)
            {
                Cancel.Register(() => input.StopRecording());
                Cancel.Register(() => tcs.TrySetCanceled());
            }

            input.StartRecording();

            var recorded_samples = await tcs.Task.ConfigureAwait(false);

            if (last_buffer_length > 0)
                input.BufferMilliseconds = last_buffer_length;

            var max_amplitude = 1 << (input.WaveFormat.BitsPerSample - 1);
            var signal_samples = new double[recorded_samples.Count];
            for (var i = 0; i < signal_samples.Length; i++)
                signal_samples[i] = (double)recorded_samples[i] / max_amplitude;

            return new SamplesDigitalSignal(1d / samples_rate, signal_samples);
        }
    }
}
