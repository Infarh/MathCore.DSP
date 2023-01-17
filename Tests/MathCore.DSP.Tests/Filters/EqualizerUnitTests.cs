using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class EqualizerUnitTests
{
    [TestMethod]
    public void SimpleTest()
    {
        const double f0 = 2e3;
        const double fd = 10e3;
        const double dt = 1 / fd;

        const double alpha = 0.0;
        const double df    = 0.5e3;

        const double w0 = Consts.pi2 * f0 * dt;
        const double dw = Consts.pi2 * df * dt;

        var equalizer0 = new EqualizerUnit(w0, dw) { Alpha = alpha };
        var equalizer = new EqualizerUnit(dt, f0, df) { Alpha = alpha };

        //var kk     = equalizer.FrequencyResponse(0, fd / 2, 100, dt).ToArray();
        //var kk_abs = kk.ToAbs();
        //var kk_arg = kk.ToArg();

        var (k_w0_re, k_w0_im) = equalizer.FrequencyResponse(0.4 * PI / Consts.pi2);
        var (k_f0_re, k_f0_im) = equalizer.FrequencyResponse(f0, dt);

        k_w0_re.AssertEquals(alpha, 2.23e-16);
        k_w0_im.AssertEquals(0);

        k_f0_re.AssertEquals(alpha, 2.23e-16);
        k_f0_im.AssertEquals(0);

        var k_fmin  = equalizer.FrequencyResponse(equalizer.fmin, dt).Abs;

        var k_fmax2 = equalizer.FrequencyResponse(equalizer.fmax + (equalizer.fmax - f0), dt).Abs;

        var kk  = alpha / k_fmin;

        var k_min_max_expected = (1 - alpha).Abs() / Consts.sqrt_2;
    }
}
