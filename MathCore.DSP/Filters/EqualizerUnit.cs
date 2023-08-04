namespace MathCore.DSP.Filters;

public class EqualizerUnit(double w0, double dw) : AllPassFilter2(-Math.Cos(w0), (1 - Math.Sin(dw)) / Math.Cos(dw))
{
    public double Alpha { get; set; } = 1;

    public double w0 { get; } = w0;

    public double dw { get; } = dw;

    public double dt { get; } = 1;

    public double f0 => w0 / (Consts.pi2 * dt);

    public double df => dw / (Consts.pi2 * dt);

    public double fmin => (w0 - 0.5 * dw) / (Consts.pi2 * dt);

    public double fmax => (w0 + 0.5 * dw) / (Consts.pi2 * dt);

    public EqualizerUnit(double dt, double f0, double df) : this(Consts.pi2 * f0 * dt, Consts.pi2 * df * dt) => this.dt = dt;

    public override Complex FrequencyResponse(double f)
    {
        var H = base.FrequencyResponse(f);
        var alpha = Alpha * 0.5;
        return (0.5 + alpha + (0.5 - alpha) * H);
    }

    public override double Process(double Sample, double[] state)
    {
        var v = base.Process(Sample, state);
        return 0.5 * (Sample + v + (Sample - v) * Alpha);
    }
}
