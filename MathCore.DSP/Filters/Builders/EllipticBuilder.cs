namespace MathCore.DSP.Filters.Builders;

public readonly ref struct EllipticBuilder
{
    public FrequencyPassType PassType { get; init; }

    public double dt { get; init; }

    public double fd { get => 1 / dt; init => dt = 1 / value; }

    public EllipticBuilder LowPass => this with { PassType = FrequencyPassType.LowPass };
    public EllipticBuilder HighPass => this with { PassType = FrequencyPassType.HighPass };
    public EllipticBuilder BandPass => this with { PassType = FrequencyPassType.BandPass };
    public EllipticBuilder BandStop => this with { PassType = FrequencyPassType.BandStop };

    private readonly double? _PassFrequency;
    private readonly double? _StopFrequency;
    private readonly double? _PassHighFrequency;
    private readonly double? _StopHighFrequency;

    private readonly int? _Order;

    private readonly double? _Gp;
    private readonly double? _Gs;
    private readonly double? _Rp;
    private readonly double? _Rs;

    public double? PassFrequency
    {
        get => _PassFrequency;
        init => _PassFrequency = value;
    }

    public double? StopFrequency
    {
        get => _StopFrequency;
        init => _StopFrequency = value;
    }

    public double? PassHighFrequency
    {
        get => _PassHighFrequency;
        init => _PassHighFrequency = value;
    }

    public double? StopHighFrequency
    {
        get => _StopHighFrequency;
        init => _StopHighFrequency = value;
    }

    public int? Order
    {
        get => _Order;
        init => _Order = value;
    }

    public double? Gp
    {
        get => _Gp;
        init => _Gp = value;
    }

    public double? Gs
    {
        get => _Gs;
        init => _Gs = value;
    }

    public double? Rp
    {
        get => _Rp;
        init => _Rp = value;
    }

    public double? Rs
    {
        get => _Rs;
        init => _Rs = value;
    }
}