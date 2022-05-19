using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

namespace MathCore.DSP.Extensions;

public static class DigitalFilterEx
{
    public static IEnumerable<(double f, Complex H)> EnumTransmissionCoefficients(
        this DigitalFilter filter,
        double Fmin,
        double Fmax,
        double dF)
    {
        var f = Fmin;
        while (f <= Fmax)
        {
            var H = filter.GetTransmissionCoefficient(f);
            yield return (f, H);
            f += dF;
        }
    }

    public static IEnumerable<(double f, double H)> EnumTransmissionCoefficientsTransient(
        this DigitalFilter filter,
        double dt,
        double tau,
        double Fmin,
        double Fmax,
        double dF)
    {
        var skip = (int)Math.Ceiling(tau / dt);

        var f = Fmin;
        while (f <= Fmax)
        {
            var x0 = f == 0 
                ? MathEnumerableSignal.Const(dt, 10000)
                : MathEnumerableSignal.Sin(dt, f, (int)(50 / f / dt) + skip);
            var y0 = filter.ProcessIndividual(x0);

            var x = new SamplesDigitalSignal(dt, x0.Skip(skip));
            var y = new SamplesDigitalSignal(dt, y0.Skip(skip));

            var H = y.Power / (f == 0 ? 1 : x.Power);

            if (f == 400)
            {
                Console.WriteLine();
            }

            yield return (f, H);
            f += dF;
        }
    }
}
