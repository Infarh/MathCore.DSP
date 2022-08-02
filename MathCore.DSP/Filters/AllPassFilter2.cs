namespace MathCore.DSP.Filters;

public class AllPassFilter2 : IIR
{
    public AllPassFilter2(double k1, double k2)
        : base(
            B: new[] { k2, k1 * (1 + k2), 01 },
            A: new[] { 01, k1 * (1 + k2), k2 })
    {
    }
}