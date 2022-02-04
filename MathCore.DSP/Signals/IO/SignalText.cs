using System.Globalization;

namespace MathCore.DSP.Signals.IO;

public class SignalText : FileSignal
{
    public NumberFormatInfo? NumbersFormat { get; set; } = NumberFormatInfo.InvariantInfo;

    public SignalText(string FileName) : base(FileName) { }

    protected override void Initialize()
    {
        var format = NumbersFormat ?? NumberFormatInfo.InvariantInfo;

        using var reader = new StreamReader(File.OpenRead(FileName));
        _dt = double.Parse(reader.ReadLine()!, format);
        _t0 = double.Parse(reader.ReadLine()!, format);
    }

    public override IEnumerable<SignalSample> GetSamples(double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity)
    {
        var format = NumbersFormat ?? NumberFormatInfo.InvariantInfo;

        using var reader = new StreamReader(File.OpenRead(FileName));
        var (dt, t0) = (double.Parse(reader.ReadLine()!, format), double.Parse(reader.ReadLine()!, format));
        if (_dt is double.NaN)
            (_dt, _t0) = (dt, t0);

        var stream = reader.BaseStream;
        var sample_index = 0;
        while (stream.Position < stream.Length)
            if (reader.ReadLine() is { Length: > 0 } line && double.TryParse(line, NumberStyles.Any, format, out var value))
            {
                var t = t0 + sample_index * dt;
                if (t < Tmin)
                {
                    sample_index++;
                    continue;
                }

                if (t > Tmax)
                    yield break;

                yield return new(t, value);
                sample_index++;
            }
    }

    public override void SetSamples(IEnumerable<double> Samples)
    {
        if (_dt is double.NaN or <= 0)
            throw new InvalidOperationException("Не задан период дискретизации");

        var format = NumbersFormat ?? NumberFormatInfo.InvariantInfo;

        using var writer = new StreamWriter(File.Create(FileName));
        writer.WriteLine(_dt);
        writer.WriteLine(_t0);

        foreach (var sample in Samples)
            writer.WriteLine(sample.ToString(format));
    }
}