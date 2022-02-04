using System.Buffers;
using System.Drawing;

namespace MathCore.DSP.Signals.IO;

public class Signal8 : FileSignal
{
    private const int __HeaderByteLength = 8 * 4;
    private const int __SampleByteLength = 1;
    private const int __MaxValue = (1 << __SampleByteLength * 8) - 1;

    private double _Min = double.NaN;
    private double _Max = double.NaN;
    private int _SamplesCount;

    public double Min
    {
        get
        {
            if (_Min is not double.NaN)
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
            if (_Max is not double.NaN)
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
            if (_Min is not double.NaN)
                return _SamplesCount;

            Initialize();
            return _SamplesCount;
        }
    }

    public double TimeLength => SamplesCount * dt;

    public double TimeMax => t0 + TimeLength;

    public Signal8(string FileName) : base(FileName) { }

    protected override void Initialize()
    {
        using var reader = new BinaryReader(File.OpenRead(FileName));
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
        using var reader = new BinaryReader(File.OpenRead(FileName));
        var (dt, t0) = (reader.ReadDouble(), reader.ReadDouble());
        var (min, max) = (reader.ReadDouble(), reader.ReadDouble());
        if (_dt is double.NaN)
        {
            (_dt, _t0) = (dt, t0);
            (_Min, _Max) = (min, max);
        }

        var delta = max - min;
        min = _Min;
        var k = delta / __MaxValue;

        const int buffer_size = 1024;
        var buffer = ArrayPool<byte>.Shared.Rent(buffer_size);
        try
        {
            var stream = reader.BaseStream;
            int readed;
            var sample_index = 0;
            do
            {
                readed = stream.Read(buffer, 0, buffer_size);
                for (var i = 0; i < readed; i++)
                {
                    var t = t0 + sample_index * dt;
                    if (t > Tmax)
                        yield break;

                    var code = buffer[i];
                    var sample = code * k + min;

                    yield return new(t, sample);
                    sample_index++;
                }
            }
            while (readed == buffer_size);
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

        using var writer = new BinaryWriter(File.Create(FileName));
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
            var code = (byte)((sample - min) * k);
            writer.Write(code);
            count++;
        }

        _SamplesCount = count;
    }
}