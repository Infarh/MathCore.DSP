using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

using static System.Math;

using static MathCore.Consts;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class LowPassRCTests : UnitTest
{
    [TestMethod]
    public void CreationTest()
    {
        const double fd = 100;
        const double dt = 1 / fd;
        const double f0 = 10;
        //const double eps = 3e-3;

        var rc = new RCLowPass(dt, f0);

        var a = rc.A;
        var b = rc.B;

        if (a is not [var a0, var a1]) throw new AssertFailedException();
        if (b is not [var b0, var b1]) throw new AssertFailedException();

        var w0 = 1 / Tan(PI * f0 * dt);
        Assert.AreEqual(1 / (1 + w0), b0 / a0);
        Assert.AreEqual(1 / (1 + w0), b1 / a0);
        //Assert.AreEqual(1, a0 / a0);
        Assert.AreEqual((1 - w0) / (1 + w0), a1 / a0);
    }

    [TestMethod]
    public void ImpulseResponseTest()
    {
        const double fd = 100;
        const double dt = 1 / fd;
        const double f0 = 10;

        var rc = new RCLowPass(dt, f0);

        var a = rc.A;
        var b = rc.B;

        var a0 = a[0];
        var a1 = a[1] / a0;
        var b0 = b[0] / a0;
        var b1 = b[1] / a0;

        List<double> expected_impulse_response =
        [
            b0,
            b0 * -a1 + b1,
            b0 * -a1 * -a1 + b1 * -a1,
            b0 * -a1 * -a1 * -a1 + b1 * -a1 * -a1,
            b0 * -a1 * -a1 * -a1 * -a1 + b1 * -a1 * -a1 * -a1,
            b0 * -a1 * -a1 * -a1 * -a1 * -a1 + b1 * -a1 * -a1 * -a1 * -a1,
        ];

        var delta = new double[expected_impulse_response.Count];
        delta[0] = 1;

        var impulse_response = rc.Process(delta).ToArray();

        const double eps = 1e-17;
        CollectionAssert.AreEqual(expected_impulse_response, impulse_response, GetComparer(eps));
    }

    [TestMethod]
    public void ProcessSignalAtCutOffFrequencyTestTest()
    {
        const double fd = 100;
        const double dt = 1 / fd;
        const double f0 = 10;
        const double eps = 3e-3;

        const double A0 = sqrt_2;
        var x0 = new SamplesDigitalSignal(dt, 1024, t => A0 * Cos(2 * PI * f0 * t));

        var x0_power = x0.Power;
        Assert.AreEqual(1, x0_power, eps);

        var rc = new RCLowPass(dt, f0);

        var y0 = rc.Process(x0);

        var y0_power = y0.Power;

        Assert.AreEqual(x0_power / 2, y0_power, eps);
    }

    [TestMethod]
    public void TransmissionCoefficientAtZeroFrequency()
    {
        const double fd = 100;
        const double dt = 1 / fd;
        const double f0 = 10;
        //const double eps = 3.02e-4;
        const double f = 0;

        var rc = new RCLowPass(dt, f0);

        var c = rc.FrequencyResponse(f, dt);

        Assert.That.Value(c.Abs).IsEqual(1);
    }

    [TestMethod]
    public void TransmissionCoefficientAtCutoffFrequency()
    {
        const double fd = 100;
        const double dt = 1 / fd;
        const double f0 = 10;
        const double eps = 1e-14;
        const double f = f0;

        var rc = new RCLowPass(dt, f0);

        var c = rc.FrequencyResponse(f, dt);

        Assert.That.Value(c.Abs).IsEqual(sqrt_2_inv, eps);
        Assert.That.Value(c.Arg * ToDeg).IsEqual(-45, eps);
    }

    [TestMethod]
    public void SpectrumTest()
    {
        const double fd = 100;
        const double dt = 1 / fd;
        const double f0 = 10;
        //const double eps = 1e-14;
        //const double f = f0;
        const int samples_count = 100;
        //const double total_time = dt * samples_count;
        //const double df = 1 / total_time;

        var rc = new RCLowPass(dt, f0);

        const double A1 = 1;
        const double f1 = 3;
        const double phi1 = 30 * ToRad;
        var s1 = new SamplesDigitalSignal(dt, samples_count, t => A1 * Cos(2 * PI * f1 * t + phi1));

        const double A2 = 8;
        const double f2 = 10;
        const double phi2 = 0 * ToRad;
        var s2 = new SamplesDigitalSignal(dt, samples_count, t => A2 * Cos(2 * PI * f2 * t + phi2));

        const double A3 = 4;
        const double f3 = 15;
        const double phi3 = 60 * ToRad;
        var s3 = new SamplesDigitalSignal(dt, samples_count, t => A3 * Cos(2 * PI * f3 * t + phi3));

        //var s = s1 + s2 + s3;
    }
}