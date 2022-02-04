namespace MathCore.DSP.Signals.IO;

public class SignalDouble : FileSignal
{
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

    public SignalDouble(string FileName) : base(FileName) { }

    protected override void Initialize()
    {
        using var reader = new BinaryReader(File.OpenRead(FileName));
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
        using var reader = new BinaryReader(File.OpenRead(FileName));
        var (dt, t0) = (reader.ReadDouble(), reader.ReadDouble());
        if (_dt is double.NaN)
            (_dt, _t0) = (dt, t0);

        var stream = reader.BaseStream;
        var sample_index = 0;
        while (stream.Position < stream.Length)
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

        using var writer = new BinaryWriter(File.Create(FileName));
        writer.Write(_dt);
        writer.Write(_t0);

        foreach (var sample in Samples)
            writer.Write(sample);
    }
}