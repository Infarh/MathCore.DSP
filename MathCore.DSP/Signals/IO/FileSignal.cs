using System.Buffers;
using System.Collections;

using MathCore.DSP.Infrastructure;

namespace MathCore.DSP.Signals.IO;

public abstract class FileSignal : IEnumerable<double>, ICollection
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
        if (dt <= 0)
            throw new ArgumentOutOfRangeException(nameof(dt), dt, "Период дискретизации должен быть больше 0");

        if (Tmin <= t0) return true;

        if (DataStream.CanSeek)
        {
            // Если поток поддерживает перемотку, то можно определить его длину
            // и можно определить доступное число отсчётов для заданного минимального времени

            var samples_count = (DataStream.Length - __HeaderByteLength) / BytePerSample;
            if (Tmin >= samples_count * dt + t0)
                return false;

            var bytes_offset = __HeaderByteLength + (long)((Tmin - t0) / dt) * BytePerSample;
            if (DataStream.Length - bytes_offset < BytePerSample)
                return false;

            DataStream.Seek(bytes_offset, SeekOrigin.Begin);
        }
        else
        {
            const int samples_count = 512;
            var buffer_size = samples_count * BytePerSample;
            var bytes_offset = (long)((Tmin - t0) / dt) * BytePerSample;
            var buffer = ArrayPool<byte>.Shared.Rent(buffer_size);
            try
            {
                while (bytes_offset > 0)
                {
                    var readed = DataStream.Read(buffer, 0, (int)Math.Min(bytes_offset, buffer_size));
                    if(readed == 0) break;
                    bytes_offset -= readed;
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
    protected readonly FileInfo _File;

    public string FilePath => _File.FullName;
    public string FileName => _File.Name;

    public bool FileExists
    {
        get
        {
            _File.Refresh();
            return _File.Exists;
        }
    }

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

    protected FileSignal(FileInfo File) => _File = File;

    protected FileSignal(string FilePath) : this(new FileInfo(FilePath)) { }

    protected abstract void Initialize();

    public abstract IEnumerable<SignalSample> GetSamples(double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity);

    public abstract void SetSamples(IEnumerable<double> Samples);

    public virtual void SetSignal(DigitalSignal Signal)
    {
        dt = Signal.dt;
        t0 = Signal.t0;
        SetSamples(Signal);
    }

    public abstract IEnumerable<double> GetValues(double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity);

    public DigitalSignal GetSignal(double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity) => 
        new SamplesDigitalSignal(dt, GetValues(Tmin, Tmax), Math.Max(Tmin, t0));

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    public IEnumerator<double> GetEnumerator() => GetSamples().Select(s => s.Value).GetEnumerator();

    public virtual int GetSamplesCount(double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity) => GetSamples(Tmin, Tmax).Count();

    #region ICollection

    int ICollection.Count => GetSamplesCount();

    bool ICollection.IsSynchronized => false;

    object? ICollection.SyncRoot => null;

    void ICollection.CopyTo(Array array, int index) => ((ICollection)GetSignal()).CopyTo(array, index);

    #endregion
}