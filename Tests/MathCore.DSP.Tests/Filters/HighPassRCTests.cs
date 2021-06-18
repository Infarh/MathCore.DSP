using System.Collections.Generic;
using System.Linq;
using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using static System.Math;

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class HighPassRCTests : UnitTest
    {
        [TestMethod]
        public void CreationTest()
        {
            const double fd = 100;
            const double dt = 1 / fd;
            const double f0 = 10;

            var rc = new HighPassRC(f0, dt);

            var a = rc.A;
            var b = rc.B;

            Assert.AreEqual(2, a.Count);
            Assert.AreEqual(2, b.Count);

            var a0 = a[0];
            var a1 = a[1];
            var b0 = b[0];
            var b1 = b[1];

            var w0 = Tan(PI * f0 * dt);
            Assert.AreEqual(1 / (w0 + 1), b0 / a0);
            Assert.AreEqual(-1 / (w0 + 1), b1 / a0);
            //Assert.AreEqual(1, a0 / a0);
            Assert.AreEqual((w0 - 1) / (w0 + 1), a1 / a0);

            Assert.That.Value(b0).IsEqual(1);
            Assert.That.Value(b1).IsEqual(-1);
            Assert.That.Value(a0).IsEqual(w0 + 1);
            Assert.That.Value(a1).IsEqual(w0 - 1);
        }

        [TestMethod]
        public void ImpulseResponseTest()
        {
            const double fd = 100;
            const double dt = 1 / fd;
            const double f0 = 10;

            var rc = new HighPassRC(f0, dt);

            var a = rc.A;
            var b = rc.B;

            var a0 = a[0];
            var a1 = a[1] / a0;
            var b0 = b[0] / a0;
            var b1 = b[1] / a0;

            var expected_impulse_response = new List<double>
            {
                b0,
                b0 * -a1 + b1,
                b0 * -a1 * -a1 + b1 * -a1,
                b0 * -a1 * -a1 * -a1 + b1 * -a1 * -a1,
                b0 * -a1 * -a1 * -a1 * -a1 + b1 * -a1 * -a1 * -a1,
                b0 * -a1 * -a1 * -a1 * -a1 * -a1 + b1 * -a1 * -a1 * -a1 * -a1,
            };

            var delta = new double[expected_impulse_response.Count];
            delta[0] = 1;

            var impulse_response = rc.Process(delta).ToArray();

            const double eps = 1.11023e-16;
            CollectionAssert.AreEqual(expected_impulse_response, impulse_response, GetComparer(eps));
        }

        [TestMethod]
        public void ProcessSignalAtCutOffFrequencyTestTest()
        {
            const double fd = 100;
            const double dt = 1 / fd;
            const double f0 = 10;
            const double eps = 3.02e-4;

            const double A0 = Consts.sqrt_2;
            var x0 = new SamplesDigitalSignal(dt, 1024, t => A0 * Cos(2 * PI * f0 * t));

            var x0_power = x0.Power;
            Assert.AreEqual(1, x0_power, eps);

            var rc = new HighPassRC(f0, dt);

            var y0 = rc.Process(x0);

            var y0_power = y0.Power;

            Assert.AreEqual(x0_power / 2, y0_power, eps);
        }

        [TestMethod]
        public void TransmissionCoefficientAtZerroFrequency()
        {
            const double fd = 100;
            const double dt = 1 / fd;
            const double f0 = 10;
            //const double eps = 3.02e-4;
            const double f = 0;

            var rc = new HighPassRC(f0, dt);

            var c = rc.GetTransmissionCoefficient(f, dt);

            Assert.That.Value(c.Abs).IsEqual(0);
        }

        [TestMethod]
        public void TransmissionCoefficientAtCutoffFrequency()
        {
            const double fd = 100;
            const double dt = 1 / fd;
            const double f0 = 10;
            //const double eps = 3.02e-4;
            const double f = f0;

            var rc = new HighPassRC(f0, dt);

            var c = rc.GetTransmissionCoefficient(f, dt);

            Assert.That.Value(c.Abs).IsEqual(Consts.sqrt_2_inv, 1e-14);
            Assert.That.Value(c.Arg * Consts.ToDeg).IsEqual(45);
        }
    }
}