using System.Buffers;
using System.Globalization;

using MathCore.Values;

namespace MathCore.DSP.Signals.IO;

public class SignalText : FileSignal
{
    public NumberFormatInfo? NumbersFormat { get; set; } = NumberFormatInfo.InvariantInfo;

    public SignalText(FileInfo File) : base(File) { }

    public SignalText(string FilePath) : base(FilePath) { }

    protected override void Initialize()
    {
        var format = NumbersFormat ?? NumberFormatInfo.InvariantInfo;

        using var reader = _File.OpenText();
        _dt = double.Parse(reader.ReadLine()!, format);
        _t0 = double.Parse(reader.ReadLine()!, format);
    }

    public override IEnumerable<SignalSample> GetSamples(double Tmin = double.NegativeInfinity, double Tmax = double.PositiveInfinity)
    {
        if (Tmax < Tmin)
            yield break;

        using var stream = _File.OpenRead();
        foreach (var sample in GetSamples(stream, Tmin, Tmax))
            yield return sample;
    }

    public static IEnumerable<SignalSample> GetSamples(Stream DataStream,
        double Tmin = double.NegativeInfinity,
        double Tmax = double.PositiveInfinity,
        NumberFormatInfo? NumbersFormat = null)
    {
        var format = NumbersFormat ?? NumberFormatInfo.InvariantInfo;

        using var reader = new StreamReader(DataStream);
        var (dt, t0) = (double.Parse(reader.ReadLine()!, format), double.Parse(reader.ReadLine()!, format));

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

        using var writer = _File.CreateText();
        writer.WriteLine(_dt);
        writer.WriteLine(_t0);

        foreach (var sample in Samples)
            writer.WriteLine(sample.ToString(format));
    }

    public override IEnumerable<double> GetValues(double Tmin = Double.NegativeInfinity, double Tmax = Double.PositiveInfinity)
    {
        if (Tmax < Tmin)
            yield break;

        using var stream = _File.OpenRead();
        foreach (var sample in GetValues(stream, Tmin, Tmax))
            yield return sample;
    }

    public static IEnumerable<double> GetValues(
        Stream DataStream, 
        double Tmin = double.NegativeInfinity,
        double Tmax = double.PositiveInfinity,
        NumberFormatInfo? NumbersFormat = null)
    {
        var format = NumbersFormat ?? NumberFormatInfo.InvariantInfo;

        using var reader = new StreamReader(DataStream);
        var (dt, t0) = (double.Parse(reader.ReadLine()!, format), double.Parse(reader.ReadLine()!, format));

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

                yield return value;
                sample_index++;
            }
    }
}