using System.Net.Security;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class EqualizerUnitTests
{
    [TestMethod]
    public void FrequencyResponseTest()
    {
        const double f0 = 2e3;
        const double fd = 10e3;
        const double dt = 1 / fd;

        const double alpha = 0.0;
        const double df    = 0.5e3;

        //const double w0 = Consts.pi2 * f0 * dt;
        //const double dw = Consts.pi2 * df * dt;

        //var equalizer0 = new EqualizerUnit(w0, dw) { Alpha = alpha };
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

    [DataTestMethod]
    [DataRow(2.0, 1.12e-14, DisplayName = "SignalProcessing_f0 alpha=2.0 eps=1.12e-14")]
    [DataRow(1.9, 1.05e-14, DisplayName = "SignalProcessing_f0 alpha=1.9 eps=1.05e-14")]
    [DataRow(1.8, 9.56e-15, DisplayName = "SignalProcessing_f0 alpha=1.8 eps=9.56e-15")]
    [DataRow(1.7, 8.67e-15, DisplayName = "SignalProcessing_f0 alpha=1.7 eps=8.67e-15")]
    [DataRow(1.6, 7.34e-15, DisplayName = "SignalProcessing_f0 alpha=1.6 eps=7.34e-15")]
    [DataRow(1.5, 6.45e-15, DisplayName = "SignalProcessing_f0 alpha=1.5 eps=6.45e-15")]
    [DataRow(1.4, 5.34e-15, DisplayName = "SignalProcessing_f0 alpha=1.4 eps=5.34e-15")]
    [DataRow(1.3, 4.23e-15, DisplayName = "SignalProcessing_f0 alpha=1.3 eps=4.23e-15")]
    [DataRow(1.2, 3.12e-15, DisplayName = "SignalProcessing_f0 alpha=1.2 eps=3.12e-15")]
    [DataRow(1.1, 1.79e-15, DisplayName = "SignalProcessing_f0 alpha=1.1 eps=1.79e-15")]
    [DataRow(1.0, 2.23e-16, DisplayName = "SignalProcessing_f0 alpha=1.0 eps=2.23e-16")]
    [DataRow(0.9, 2.12e-15, DisplayName = "SignalProcessing_f0 alpha=0.9 eps=2.12e-15")]
    [DataRow(0.8, 4.01e-15, DisplayName = "SignalProcessing_f0 alpha=0.8 eps=4.01e-15")]
    [DataRow(0.7, 6.12e-15, DisplayName = "SignalProcessing_f0 alpha=0.7 eps=6.12e-15")]
    [DataRow(0.6, 7.89e-15, DisplayName = "SignalProcessing_f0 alpha=0.6 eps=7.89e-15")]
    [DataRow(0.5, 1.02e-14, DisplayName = "SignalProcessing_f0 alpha=0.5 eps=1.02e-14")]
    [DataRow(0.4, 1.24e-14, DisplayName = "SignalProcessing_f0 alpha=0.4 eps=1.24e-14")]
    [DataRow(0.3, 1.46e-14, DisplayName = "SignalProcessing_f0 alpha=0.3 eps=1.46e-14")]
    [DataRow(0.2, 1.68e-14, DisplayName = "SignalProcessing_f0 alpha=0.2 eps=1.68e-14")]
    [DataRow(0.1, 1.90e-14, DisplayName = "SignalProcessing_f0 alpha=0.1 eps=1.90e-14")]
    [DataRow(0.0, 5.17e-14, DisplayName = "SignalProcessing_f0 alpha=0.0 eps=5.17e-14")]
    public void SignalProcessing_f0(double alpha, double eps)
    {
        const double f0 = 2e3;
        const double fd = 10e3;
        const double dt = 1 / fd;

        const double df = 0.5e3;

        var equalizer = new EqualizerUnit(dt, f0, df) { Alpha = alpha };

        const int samples_count_period_f0 = (int)(fd / f0);
        const int samples_count = samples_count_period_f0 * 100;
        var x_f0 = Signal.Harmonic.Sin(f0, dt, samples_count);
        var x_f0_p2p = x_f0.PeakToPeakAmplitude;

        var y_f0 = equalizer.ProcessIndividual(x_f0);
        var y1_f0 = new EnumerableSignal(dt, y_f0.Skip(samples_count_period_f0 * 50));

        var y_f0_p2p = y1_f0.PeakToPeakAmplitude;

        var k_f0_p2p = y_f0_p2p / x_f0_p2p;

        k_f0_p2p.AssertEquals(alpha, eps);
   }
}
