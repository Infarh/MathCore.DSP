using System.Buffers;
using System;
using System.Collections;

using MathCore.DSP.Infrastructure;

namespace MathCore.DSP.Signals.IO;

public abstract class FileSignal : IEnumerable<double>
{
    protected const int __HeaderByteLength = 8 * 4;

    protected static (double dt, double t0, double Min, double Max) ReadHeader(Stream DataStream)
    {
        Span<byte> buffer = stackalloc byte[__HeaderByteLength];
        DataStream.ReadToSpan(buffer);
        var values = buffer.Cast<double>();
        return (values[0], values[1], values[2], values[3]);
    }

    protected static bool FindTmin(Stream DataStream, double Tmin, double dt, double t0, int BytePerSample)
    {
        if (Tmin <= t0) return true;

        if (DataStream.CanSeek)
        {
            // Если поток поддерживает перемотку, то можно определить его длину
            // и можно определить доступное число отсчётов для заданного минимального времени

            var samples_count = (DataStream.Length - __HeaderByteLength) / BytePerSample;
            if (Tmin >= samples_count * dt + t0)
                return false;

            var bytes_offset = (long)((Tmin - t0) / dt) + __HeaderByteLength;
            if (bytes_offset > DataStream.Length)
                return false;

            DataStream.Seek(bytes_offset, SeekOrigin.Begin);
        }
        else
        {
            const int samples_count = 512;
            var buffer_size = samples_count * BytePerSample;
            var bytes_offset = (int)((Tmin - t0) / BytePerSample);
            var buffer = ArrayPool<byte>.Shared.Rent(buffer_size);
            try
            {
                while (bytes_offset >= buffer_size)
                {
                    _ = DataStream.Read(buffer, 0, Math.Min(bytes_offset, buffer_size));
                    bytes_offset -= buffer_size;
                }
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(buffer);
            }
        }

        return true;
    }

    protected double _dt = double.NaN;
    protected double _t0;

    public string FileName { get; }

    public double dt
    {
        get
        {
            if (_dt > 0) return _dt;
            Initialize();
            return _dt;
        }
        set => _dt = value;
    }

    public double t0
    {
        get
        {
            if (_dt > 0) return _t0;
            Initialize();
            return _t0;
        }
        set => _t0 = value;
    }

    protected FileSignal(string FileName) => this.FileName = FileName;

    protected abstract void Initialize();

    public abstract IEnumerable<SignalSample> GetSamples(double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity);

    public abstract void SetSamples(IEnumerable<double> Samples);

    public DigitalSignal GetSignal() => new SamplesDigitalSignal(dt, GetSamples().Select(s => s.Value));

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    public IEnumerator<double> GetEnumerator() => GetSamples().Select(s => s.Value).GetEnumerator();
}