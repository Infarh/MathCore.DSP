using System.Collections.Generic;
using System.Linq;
using MathCore.DSP.Filters;
using MathCore.DSP.Fourier;
using MathCore.DSP.Signals;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using static System.Math;

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class BandPassRLCTests : UnitTest
    {
        [TestMethod]
        public void CreationTest()
        {
            const double fd = 1000;
            const double dt = 1 / fd;
            const double f0 = 30;
            const double Df = 10;

            var w = Tan(PI * f0 * dt);
            var dw = PI * Df * dt;

            var expected_b0 = dw;
            var expected_b1 = 0d;
            var expected_b2 = -dw;

            var expected_a0 = w * w + dw + 1;
            var expected_a1 = 2 * (w * w - 1);
            var expected_a2 = w * w - dw + 1;

            var rlc = new BandPassRLC(f0, Df, dt);

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
            const double Df = 10;

            var rlc = new BandPassRLC(f0, Df, dt);

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

            const double eps = 3e-17;
            CollectionAssert.AreEqual(expected_impulse_response, impulse_response, GetComparer(eps));
        }

        [TestMethod]
        public void ProcessSignalAtCentralFrequencyTestTest()
        {
            const double fd = 1000;
            const double dt = 1 / fd;
            const double f0 = 30;
            const double Df = 10;
            const double eps = 0.048;

            const double A0 = Consts.sqrt_2;
            var x0 = new SamplesDigitalSignal(dt, 1024, t => A0 * Cos(2 * PI * f0 * t));

            var x0_power = x0.Power;
            Assert.AreEqual(1, x0_power, eps);

            var rlc = new BandPassRLC(f0, Df, dt);

            var y0 = rlc.Process(x0);

            var y0_power = y0.Power;

            Assert.AreEqual(x0_power, y0_power, eps, $"delta:{Abs(x0_power - y0_power)}");
        }

        [TestMethod]
        public void ProcessSignalAtZerroFrequencyTestTest()
        {
            const double fd = 1000;
            const double dt = 1 / fd;
            const double f0 = 30;
            const double Df = 10;
            const double eps = 0.86e-3;

            var x0 = new SamplesDigitalSignal(dt, 1024, t => 1 * Cos(2 * PI * 0 * t));

            var x0_power = x0.Power;
            Assert.AreEqual(1, x0_power, eps);

            var rlc = new BandPassRLC(f0, Df, dt);

            var y0 = rlc.Process(x0);

            var y0_power = y0.Power;

            Assert.AreEqual(0, y0_power, eps, $"delta:{Abs(0 - y0_power)}");
        }

        [TestMethod]
        public void ProcessSignalAtHighFrequencyTestTest()
        {
            const double fd = 1000;
            const double dt = 1 / fd;
            const double f0 = 30;
            const double Df = 10;
            const double eps = 7.7e-6;

            var x0 = new SamplesDigitalSignal(dt, 1024, t => 1 * Cos(2 * PI * fd / 2 * t));

            var x0_power = x0.Power;
            Assert.AreEqual(1, x0_power, eps);

            var rlc = new BandPassRLC(f0, Df, dt);

            var y0 = rlc.Process(x0);

            var y0_power = y0.Power;

            Assert.AreEqual(0, y0_power, eps, $"delta:{Abs(0 - y0_power):e2}");
        }

        [TestMethod]
        public void TransmissionCoefficientAtCentralFrequencyTestTest()
        {
            const double fd = 1000;
            const double dt = 1 / fd;
            const double f0 = 30;
            const double Df = 10;

            var rlc = new BandPassRLC(f0, Df, dt);

            var c = rlc.GetTransmissionCoefficient(f0, dt);

            const double eps = 1.12e-15;
            Assert.AreEqual(1, c.Abs, eps, $"delta:{Abs(1 - c.Abs):e2}");
        }

        [TestMethod]
        public void TransmissionCoefficientAtZerroFrequencyTestTest()
        {
            const double fd = 1000;
            const double dt = 1 / fd;
            const double f0 = 30;
            const double Df = 10;

            var rlc = new BandPassRLC(f0, Df, dt);

            var c = rlc.GetTransmissionCoefficient(0, dt);

            Assert.AreEqual(0, c.Abs, $"delta:{Abs(0 - c.Abs):e2}");
        }

        [TestMethod]
        public void TransmissionCoefficientAtHighFrequencyTestTest()
        {
            const double fd = 1000;
            const double dt = 1 / fd;
            const double f0 = 30;
            const double Df = 10;

            var rlc = new BandPassRLC(f0, Df, dt);

            var c = rlc.GetTransmissionCoefficient(fd / 2, dt);

            const double eps = 1.93e-18;
            Assert.AreEqual(0, c.Abs, eps, $"delta:{Abs(0 - c.Abs):e2}");
        }

        [TestMethod]
        public void SpectrumPassTest()
        {
            const double fd = 1000;
            const double dt = 1 / fd;
            const double f0 = 30;
            const double Df = 10;

            var rlc = new BandPassRLC(f0, Df, dt);

            const int samples_count = 1000;
            var s0 = new SamplesDigitalSignal(dt, Enumerable.Repeat(1, samples_count));
            var s1 = new SamplesDigitalSignal(dt, samples_count, t => 1 * Cos(2 * PI * (f0 - Df / 2) * t));
            var s2 = new SamplesDigitalSignal(dt, samples_count, t => 1 * Cos(2 * PI * f0 * t));
            var s3 = new SamplesDigitalSignal(dt, samples_count, t => 1 * Cos(2 * PI * (f0 + Df / 2) * t));

            var s = s0 + s1 + s2 + s3;
            s.FourierTransform().ToAbsArg(out var abs_S, out var arg_S);
            var abs_S_25 = abs_S[25];
            var abs_S_30 = abs_S[30];
            var abs_S_35 = abs_S[35];

            var arg_S_25 = arg_S[25];
            var arg_S_30 = arg_S[30];
            var arg_S_35 = arg_S[35];

            var y = rlc.ProcessIndividual(s);
            y.FourierTransform().ToAbsArg(out var abs_Y, out var arg_Y);
            Assert.That.Value(abs_Y[0]).IsEqual(0, 1.68e-3);
            var abs_Y_25 = abs_Y[25];
            var abs_Y_30 = abs_Y[30];
            var abs_Y_35 = abs_Y[35];

            var H_25 = rlc.GetTransmissionCoefficient(f0 - Df / 2, dt);
            var H_30 = rlc.GetTransmissionCoefficient(f0, dt);
            var H_35 = rlc.GetTransmissionCoefficient(f0 + Df / 2, dt);

            var H_25_abs = H_25.Abs;
            var H_30_abs = H_30.Abs;
            var H_35_abs = H_35.Abs;

            var K25 = abs_Y_25 / abs_S_25;
            var K30 = abs_Y_30 / abs_S_30;
            var K35 = abs_Y_35 / abs_S_35;

            Assert.That.Value(K25).IsEqual(H_25_abs, 0.042);
            Assert.That.Value(K30).IsEqual(H_30_abs, 0.062);
            Assert.That.Value(K35).IsEqual(H_35_abs, 0.045);

            var arg_Y_25 = arg_Y[25];
            var arg_Y_30 = arg_Y[30];
            var arg_Y_35 = arg_Y[35];

            var H_25_arg = H_25.Arg;
            var H_30_arg = H_30.Arg;
            var H_35_arg = H_35.Arg;

            var Arg25 = arg_S_25 - arg_Y_25;
            var Arg30 = arg_S_30 - arg_Y_30;
            var Arg35 = arg_S_35 - arg_Y_35;

            Assert.That.Value(Arg25).IsEqual(H_25_arg, 1.82e-2);
            Assert.That.Value(Arg30).IsEqual(H_30_arg, 1.71e-2);
            Assert.That.Value(Arg35).IsEqual(H_35_arg, 1.66e-2);
        }
    }
}