using System.Collections.Generic;
using System.Linq;
using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using static System.Math;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class BandStopRLCTests : UnitTest
{
    [TestMethod]
    public void CreationTest()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const double f0 = 30;
        const double delta_f = 10;

        var w = Tan(PI * f0 * dt);
        const double dw = PI * delta_f * dt;

        var expected_b0 = w * w + 1;
        var expected_b1 = 2 * (w * w - 1);
        var expected_b2 = w * w + 1;

        var expected_a0 = w * w + dw + 1;
        var expected_a1 = 2 * (w * w - 1);
        var expected_a2 = w * w - dw + 1;

        var rlc = new BandStopRLC(f0, delta_f, dt);

        var a = rlc.A;
        var b = rlc.B;

        Assert.AreEqual(expected_b0, b[0]);
        Assert.AreEqual(expected_b1, b[1]);
        Assert.AreEqual(expected_b2, b[2]);

        Assert.AreEqual(expected_a0, a[0]);
        Assert.AreEqual(expected_a1, a[1]);
        Assert.AreEqual(expected_a2, a[2]);
    }

    [TestMethod]
    public void ImpulseResponseTest()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const double f0 = 30;
        const double delta_f = 10;

        var rlc = new BandStopRLC(f0, delta_f, dt);

        var a = rlc.A;
        var b = rlc.B;

        var a0 = a[0];
        var a1 = a[1] / a0;
        var a2 = a[2] / a0;
        var b0 = b[0] / a0;
        var b1 = b[1] / a0;
        var b2 = b[2] / a0;

        var expected_impulse_response = new List<double>
        {
            b0,
            b0 * -a1 + b1,
            b0 * (-a1 * -a1 + -a2) + b1 * -a1 + b2,
            b0 * (-a1 * -a2 + -a1 * (-a1 * -a1 + -a2)) + b1 * (-a1 * -a1 + -a2) + b2 * -a1
        };

        var delta = new double[expected_impulse_response.Count];
        delta[0] = 1;

        var impulse_response = rlc.Process(delta).ToArray();

        const double eps = 5.6e-16;
        CollectionAssert.AreEqual(expected_impulse_response, impulse_response, GetComparer(eps));
    }

    [TestMethod]
    public void ProcessSignalAtCentralFrequencyTestTest()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const double f0 = 30;
        const double delta_f = 10;
        const double eps = 0.048;

        const double A0 = Consts.sqrt_2;
        var x0 = new SamplesDigitalSignal(dt, 1024, t => A0 * Cos(2 * PI * f0 * t));

        var x0_power = x0.Power;
        Assert.AreEqual(1, x0_power, eps);

        var rlc = new BandStopRLC(f0, delta_f, dt);

        var y0 = rlc.Process(x0);

        var y0_power = y0.Power;

        Assert.AreEqual(0, y0_power, eps, $"delta:{Abs(0 - y0_power)}");
    }

    [TestMethod]
    public void ProcessSignalAtZerroFrequencyTestTest()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const double f0 = 30;
        const double delta_f = 10;
        const double eps = 2.59e-3;

        var x0 = new SamplesDigitalSignal(dt, 1024, t => 1 * Cos(2 * PI * 0 * t));

        var x0_power = x0.Power;
        Assert.AreEqual(1, x0_power, eps);

        var rlc = new BandStopRLC(f0, delta_f, dt);

        var y0 = rlc.Process(x0);

        var y0_power = y0.Power;

        Assert.AreEqual(x0_power, y0_power, eps, $"delta:{Abs(x0_power - y0_power):e2}");
    }

    [TestMethod]
    public void ProcessSignalAtHighFrequencyTestTest()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const double f0 = 30;
        const double delta_f = 10;
        const double eps = 2.31e-5;

        var x0 = new SamplesDigitalSignal(dt, 1024, t => 1 * Cos(2 * PI * fd / 2 * t));

        var x0_power = x0.Power;
        Assert.AreEqual(1, x0_power, eps);

        var rlc = new BandStopRLC(f0, delta_f, dt);

        var y0 = rlc.Process(x0);

        var y0_power = y0.Power;

        Assert.AreEqual(x0_power, y0_power, eps, $"delta:{Abs(x0_power - y0_power):e2}");
    }

    [TestMethod]
    public void TransmissionCoefficientAtCentralFrequencyTestTest()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const double f0 = 30;
        const double delta_f = 10;

        var rlc = new BandStopRLC(f0, delta_f, dt);

        var c = rlc.FrequencyResponse(f0, dt);

        const double eps = 1.91e-14;
        Assert.AreEqual(0, c.Abs, eps, $"delta:{Abs(0 - c.Abs):e2}");
    }

    [TestMethod]
    public void TransmissionCoefficientAtZeroFrequencyTestTest()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const double f0 = 30;
        const double delta_f = 10;

        var rlc = new BandStopRLC(f0, delta_f, dt);

        var c = rlc.FrequencyResponse(0, dt);

        Assert.AreEqual(1, c.Abs, $"delta:{Abs(1 - c.Abs):e2}");
    }

    [TestMethod]
    public void TransmissionCoefficientAtHighFrequencyTestTest()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const double f0 = 30;
        const double delta_f = 10;

        var rlc = new BandStopRLC(f0, delta_f, dt);

        var c = rlc.FrequencyResponse(fd / 2, dt);

        const double eps = 1.93e-18;
        Assert.AreEqual(1, c.Abs, eps, $"delta:{Abs(1 - c.Abs):e2}");
    }
}