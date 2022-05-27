using System;
using System.Linq;
using MathCore.DSP.Filters;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class IIRTests : UnitTest
{
    [TestMethod]
    public void ImpulseResponseTest()
    {
        //             a0            0.3
        // H(z) = ------------ = ------------
        //        b0 + b1 z^-1   1 - 0.7 z^-1

        var b0 = 0.3;
        var a0 = 1;
        var a1 = -0.7;

        var b = new[] { b0 };
        var a = new[] { a0, a1 };
        var filter = new IIR(b, a);

        var impulse_response_length = 3;
        var delta = new double[3];
        delta[0] = 1;

        var impulse_response = filter.Process(delta).ToArray();

        Assert.AreEqual(impulse_response_length, impulse_response.Length);

        double[] expected_impulse_response = { b0, -a1 * b0, -a1 * -a1 * b0 };
        CollectionAssert.AreEqual(expected_impulse_response, impulse_response);
    }

    [TestMethod]
    public void SerialConnectionTest()
    {
        const double fd = 100;
        const double dt = 1 / fd;
        const double f0RC = 40;
        const double f0CR = 5;

        var RC = new LowPassRC(f0RC, dt);
        var CR = new HighPassRC(f0CR, dt);

        var filter = RC.ConnectionSerialTo(CR);

        var A = filter.A;
        var B = filter.B;

        double[] expected_a =
        {
            1.53476636079571,
            -0.33307051181674,
            -0.568158087680825
        };
        double[] expected_b = { 1, 0, -1 };
        CollectionAssert.AreEqual(expected_b, B);
        CollectionAssert.AreEqual(expected_a, A, GetComparer(1e-15));

        void CheckTransmissionCoefficient(double f)
        {
            var H = filter.FrequencyResponse(f, dt);
            var HRC = RC.FrequencyResponse(f, dt);
            var HCR = CR.FrequencyResponse(f, dt);
            var delta = HRC * HCR - H;
            const double eps = 1.25e-15;
            Assert.That.Value(delta.Abs).IsEqual(0, eps);
        }

        CheckTransmissionCoefficient(0);
        CheckTransmissionCoefficient(f0RC);
        CheckTransmissionCoefficient(f0CR);
        CheckTransmissionCoefficient(fd / 2);

        var H0 = filter.FrequencyResponse(0, dt);
        var Hf0RC = filter.FrequencyResponse(f0RC, dt);
        var Hf0CR = filter.FrequencyResponse(f0CR, dt);
        var Hfd05 = filter.FrequencyResponse(fd / 2, dt);

        Assert.AreEqual(0, H0.Abs);
        Assert.AreEqual(0, Hfd05.Abs, 1.89e-16);
        Assert.AreEqual(0, H0.Arg);
        Assert.AreEqual(-Math.PI / 2, Hfd05.Arg);

        Assert.That.Value(Hf0RC.Abs).IsEqual(Consts.sqrt_2_inv, 9.35e-4);
        Assert.That.Value(Hf0CR.Abs).IsEqual(Consts.sqrt_2_inv, 9.35e-4);
    }

    [TestMethod]
    public void ParallelConnectionTest()
    {
        const double fd = 100;
        const double dt = 1 / fd;
        const double f0RC = 40;
        const double f0CR = 5;

        var RC = new LowPassRC(f0RC, dt);
        var CR = new HighPassRC(f0CR, dt);

        var filter = RC.ConnectionParallelTo(CR);

        var A = filter.A;
        var B = filter.B;

        double[] expected_a =
        {
            1.53476636079571,
            -0.33307051181674,
            -0.568158087680825
        };
        double[] expected_b =
        {
            1.241652068278721,
            -0.16653525590837,
            -0.758347931721279
        };
        CollectionAssert.AreEqual(expected_b, B, GetComparer(1e-15));
        CollectionAssert.AreEqual(expected_a, A, GetComparer(1e-15));

        void CheckTransmissionCoefficient(double f)
        {
            var H = filter.FrequencyResponse(f, dt);
            var HRC = RC.FrequencyResponse(f, dt);
            var HCR = CR.FrequencyResponse(f, dt);
            var delta = (HRC + HCR) / 2 - H;
            const double eps = 1.25e-15;
            Assert.That.Value(delta.Abs).IsEqual(0, eps);
        }

        CheckTransmissionCoefficient(0);
        CheckTransmissionCoefficient(f0RC);
        CheckTransmissionCoefficient(f0CR);
        CheckTransmissionCoefficient(fd / 2);
    }
}