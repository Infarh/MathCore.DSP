using System.Buffers;

namespace MathCore.DSP.Signals.IO;

public class Signal8 : FileSignal
{
    private double _Min = double.NaN;
    private double _Max = double.NaN;

    public string FileName { get; }

    public double Min
    {
        get
        {
            if (_Min > 0) return _Min;
            Initialize();
            return _Min;
        }
        set => _Min = value;
    }

    public double Max
    {
        get
        {
            if (_Max > 0) return _Max;
            Initialize();
            return _Max;
        }
        set => _Max = value;
    }

    public Interval Interval { get => new(Min, Max); set => (Min, Max) = value; }

    public double Delta => _Max - _Min;

    public double Middle => (_Max + _Min) / 2;

    public Signal8(string FileName) : base(FileName) { }

    protected override void Initialize()
    {
        using var reader = new BinaryReader(File.OpenRead(FileName));
        _dt = reader.ReadDouble();
        _t0 = reader.ReadDouble();
        _Min = reader.ReadDouble();
        _Max = reader.ReadDouble();
    }

    public override IEnumerable<double> GetSamples()
    {
        using var reader = new BinaryReader(File.OpenRead(FileName));
        var (dt, t0) = (reader.ReadDouble(), reader.ReadDouble());
        var (min, max) = (reader.ReadDouble(), reader.ReadDouble());
        if (_dt is double.NaN)
        {
            (_dt, _t0) = (dt, t0);
            (_Min, _Max) = (min, max);
        }

        var delta = Delta;
        min = _Min;

        const int buffer_size = 1024;
        var buffer = ArrayPool<byte>.Shared.Rent(buffer_size);
        try
        {
            var stream = reader.BaseStream;
            int readed;
            do
            {
                readed = stream.Read(buffer, 0, buffer_size);
                for (var i = 0; i < readed; i++) 
                    yield return buffer[i] * delta + min;
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

        //foreach (var sample in Samples)
        //{
        //    var code = (byte)((sample - min) / delta);
        //    writer.Write(code);
        //}

        const int buffer_size = 1024;
        var buffer = ArrayPool<byte>.Shared.Rent(buffer_size);
        try
        {
            foreach (var chunk in Samples.AsBlockEnumerable(buffer_size))
            {
                    
            }
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(buffer);
        }
    }
}