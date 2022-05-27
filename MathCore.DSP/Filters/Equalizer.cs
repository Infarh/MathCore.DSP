namespace MathCore.DSP.Filters;

public class Equalizer : AllPassFilter2
{
    public double Alpha { get; set; } = 1;

    public Equalizer(double w0, double dw) : base(-Math.Cos(w0), (1 - Math.Sin(dw)) / Math.Cos(dw))
    {
    }

    public override Complex FrequencyResponse(double f)
    {
        var H = base.FrequencyResponse(f);
        var alpha = Alpha;
        return (1 + alpha + (1 - alpha) * H) / 2;
    }

    public override Complex FrequencyResponse(double f, double dt)
    {
        var H = base.FrequencyResponse(f, dt);
        var alpha = Alpha;
        return (1 + alpha + (1 - alpha) * H) / 2;
    }

    public override double Process(double Sample, double[] state)
    {
        var v = base.Process(Sample, state);
        return (Sample + v + (Sample - v) * Alpha) / 2;
    }
}
