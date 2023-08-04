namespace MathCore.DSP.Filters;

public static class FilterEx
{
    public static IEnumerable<Complex> FrequencyResponse(this Filter filter, double f1, double f2, double df)
    {
        var f = f1;
        while (f <= f2)
        {
            yield return filter.FrequencyResponse(f);
            f += df;
        }
    }

    public static IEnumerable<Complex> FrequencyResponse(this Filter filter, double f1, double f2, double df, double dt)
    {
        var f = f1;
        while (f <= f2)
        {
            yield return filter.FrequencyResponse(f * dt);
            f += df;
        }
    }
}
