namespace MathCore.DSP.Signals.IO;

public class SignalDouble : FileSignal
{
    public SignalDouble(string FileName) : base(FileName) { }

    protected override void Initialize()
    {
        using var reader = new BinaryReader(File.OpenRead(FileName));
        _dt = reader.ReadDouble();
        _t0 = reader.ReadDouble();
    }

    public override IEnumerable<double> GetSamples()
    {
        using var reader = new BinaryReader(File.OpenRead(FileName));
        var (dt, t0) = (reader.ReadDouble(), reader.ReadDouble());
        if (_dt is double.NaN)
            (_dt, _t0) = (dt, t0);

        var stream = reader.BaseStream;
        while (stream.Position < stream.Length)
            yield return reader.ReadDouble();
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