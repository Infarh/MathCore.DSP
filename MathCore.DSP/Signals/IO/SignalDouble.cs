namespace MathCore.DSP.Signals.IO;

public class SignalDouble : FileSignal
{
    private const int __SampleByteLength = 8;

    private int _SamplesCount;

    public int SamplesCount
    {
        get
        {
            if (_SamplesCount > 0)
                return _SamplesCount;

            Initialize();
            return _SamplesCount;
        }
    }

    public double TimeLength => SamplesCount * dt;

    public double TimeMax => t0 + TimeLength;

    public SignalDouble(FileInfo File) : base(File) { }
    public SignalDouble(string FilePath) : base(FilePath) { }

    protected override void Initialize()
    {
        using var reader = _File.OpenBinary();
        _dt = reader.ReadDouble();
        _t0 = reader.ReadDouble();

        var file_length = reader.BaseStream.Length;
        const int header_byte_length = 8 * 2;
        var data_length = file_length - header_byte_length;
        const int sample_byte_length = 8;
        _SamplesCount = (int)(data_length / sample_byte_length);
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
        using var reader = new BinaryReader(DataStream);
        var dt = reader.ReadDouble();
        var t0 = reader.ReadDouble();

        if (!FindTmin(DataStream, Tmin, dt, t0, __SampleByteLength))
            yield break;

        var sample_index = 0;
        var total_length = DataStream.Length;
        while (DataStream.Position < total_length)
        {
            var t = t0 + sample_index * dt;
            if (t > Tmax)
                yield break;

            var sample = reader.ReadDouble();
            yield return new(t, sample);
            sample_index++;
        }
    }

    public override void SetSamples(IEnumerable<double> Samples)
    {
        if (_dt is double.NaN or <= 0)
            throw new InvalidOperationException("Не задан период дискретизации");

        using var writer = _File.CreateBinary();
        writer.Write(_dt);
        writer.Write(_t0);

        foreach (var sample in Samples)
            writer.Write(sample);
    }

    public override IEnumerable<double> GetValues(double Tmin = Double.NegativeInfinity, double Tmax = Double.PositiveInfinity)
    {
        if (Tmax < Tmin)
            yield break;

        using var stream = _File.OpenRead();
        foreach (var sample in GetValues(stream, Tmin, Tmax))
            yield return sample;
    }

    public static IEnumerable<double> GetValues(Stream DataStream, double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity)
    {
        using var reader = new BinaryReader(DataStream);
        var dt = reader.ReadDouble();
        var t0 = reader.ReadDouble();

        if (!FindTmin(DataStream, Tmin, dt, t0, __SampleByteLength))
            yield break;

        var sample_index = 0;
        var total_length = DataStream.Length;
        while (DataStream.Position < total_length)
        {
            var t = t0 + sample_index * dt;
            if (t > Tmax)
                yield break;

            var sample = reader.ReadDouble();
            yield return sample;
            sample_index++;
        }
    }
}