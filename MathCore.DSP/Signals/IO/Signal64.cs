using System.Buffers;

namespace MathCore.DSP.Signals.IO;

public class Signal64 : FileSignal
{
    private const int __SampleByteLength = 8;
    private const ulong __MaxValue = ulong.MaxValue;

    private double _Min = double.NaN;
    private double _Max = double.NaN;
    private int _SamplesCount;

    public double Min
    {
        get
        {
            if (!double.IsNaN(_Min))
                return _Min;

            Initialize();
            return _Min;
        }
        set => _Min = value;
    }

    public double Max
    {
        get
        {
            if (!double.IsNaN(_Max))
                return _Max;

            Initialize();
            return _Max;
        }
        set => _Max = value;
    }

    public Interval Interval { get => new(Min, Max); set => (Min, Max) = value; }

    public double Delta => _Max - _Min;

    public double dA => Delta / __MaxValue;

    public double Middle => (_Max + _Min) / 2;

    public int SamplesCount
    {
        get
        {
            if (!double.IsNaN(_Min))
                return _SamplesCount;

            Initialize();
            return _SamplesCount;
        }
    }

    public double TimeLength => SamplesCount * dt;

    public double TimeMax => t0 + TimeLength;

    public Signal64(FileInfo File) : base(File) { }

    public Signal64(string FilePath) : base(FilePath) { }

    protected override void Initialize()
    {
        using var reader = _File.OpenBinary();
        _dt = reader.ReadDouble();
        _t0 = reader.ReadDouble();
        _Min = reader.ReadDouble();
        _Max = reader.ReadDouble();

        var file_length = reader.BaseStream.Length;
        var data_length = file_length - __HeaderByteLength;
        _SamplesCount = (int)(data_length / __SampleByteLength);
    }

    public override IEnumerable<SignalSample> GetSamples(double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity)
    {
        if (Tmax < Tmin)
            yield break;

        using var stream = _File.OpenRead();
        foreach (var sample in GetSamples(stream, Tmin, Tmax))
            yield return sample;
    }

    public static IEnumerable<SignalSample> GetSamples(Stream DataStream, double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity)
    {
        if (Tmax < Tmin) yield break;

        var (dt, t0, min, max) = ReadHeader(DataStream);

        if (!FindTmin(DataStream, Tmin, dt, t0, __SampleByteLength))
            yield break;

        var delta = max - min;
        var k = delta / __MaxValue;

        const int default_samples_count = 512;
        var samples_count = double.IsInfinity(Tmin) || double.IsInfinity(Tmax)
            ? default_samples_count
            : Math.Min(default_samples_count, (int)((Tmax - Tmin) / dt));
        var buffer_size = samples_count * __SampleByteLength;
        var buffer = ArrayPool<byte>.Shared.Rent(buffer_size);
        //var values = buffer.AsMemory().Cast<ushort>();
        try
        {
            int readed_bytes_count;
            var sample_index = 0;
            if (Tmin > t0)
            {
                t0 += (int)((Tmin - t0) / dt) * dt;
                if (t0 < Tmin)
                    t0 += dt;
            }

            do
            {
                readed_bytes_count = DataStream.Read(buffer, 0, buffer_size);
                for (var i = 0; i < readed_bytes_count; i += __SampleByteLength)
                {
                    var t = t0 + sample_index * dt;
                    if (t > Tmax)
                        yield break;

                    var code = BitConverter.ToUInt64(buffer, i);
                    var sample = code * k + min;

                    yield return new(t, sample);
                    sample_index++;
                }
            }
            while (readed_bytes_count == buffer_size);
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(buffer);
        }
    }

    public override void SetSamples(IEnumerable<double> Samples)
    {
        if (_dt is double.NaN or <= 0)
            throw new InvalidOperationException("Не задан период дискретизации");
        if (_Min is double.NaN)
            throw new InvalidOperationException("Не задан минимум интервала значений сигнала");
        if (_Max is double.NaN)
            throw new InvalidOperationException("Не задан максимум интервала значений сигнала");

        using var writer = _File.CreateBinary();
        writer.Write(_dt);
        writer.Write(_t0);

        double min, max;
        writer.Write(min = _Min);
        writer.Write(max = _Max);
        var delta = max - min;
        var k = __MaxValue / delta;

        var count = 0;
        foreach (var sample in Samples)
        {
            var code = (ulong)((sample - min) * k);
            writer.Write(code);
            count++;
        }

        _SamplesCount = count;
    }
}