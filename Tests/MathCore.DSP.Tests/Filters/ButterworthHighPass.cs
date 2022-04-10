using System;
using System.Linq;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;
using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ButterworthHighPass
{
    [TestMethod]
    public void Creation()
    {
        const double fd = 5000;                 // Гц // Частота дискретизации
        const double fs = 0.5 * fd / PI;   // Гц // Граничная частота полосы пропускания
        const double fp = 1.0 * fd / PI;   // Гц // Граничная частота полосы запирания
        const double dt = 1 / fd;               // 2с // Период дискретизации
        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 30;                   // Неравномерность в полосе пропускания (дБ)
        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        const double ws = 2 * PI * fs / fd;
        const double wp = 2 * PI * fp / fd;

        var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);
        var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);

        var Ws = 2 * PI * Fs;
        var Wp = 2 * PI * Fp;

        var eps_p = Sqrt(Pow(10, Rp / 10) - 1);
        var eps_s = Sqrt(Pow(10, Rs / 10) - 1);

        var kEps = eps_s / eps_p;
        var kW = Fp / Fs;

        kEps.AssertEquals(62.114845067566208);
        kW.AssertEquals(2.8508157176809257);

        var N = (int)Ceiling(Log(kEps) / Log(kW));
        N.AssertEquals(4);

        var L = N / 2;
        var r = N % 2;

        var alpha = eps_p.Pow(-1d / N);
        alpha.AssertEquals(1.184003988964071);

        var th0 = PI / N;

        var poles = new Complex[N];
        if (r != 0) poles[0] = -alpha;
        for (var i = r; i < poles.Length; i += 2)
        {
            var w = th0 * (i + 1 - r - 0.5);
            poles[i] = (-alpha * Sin(w), alpha * Cos(w));
            poles[i + 1] = poles[i].ComplexConjugate;
        }

        var translated_poles = poles.ToArray(p => Wp / p);
        translated_poles.AssertEquals(
            (-5033.717278134342, -12152.468522023708),
            (-5033.717278134342, +12152.468522023708),
            (-12152.468522023708, +5033.7172781343415),
            (-12152.468522023708, -5033.7172781343415)
        );


        //var z_poles = translated_poles.ToArray(p => DigitalFilter.ToZ(p, dt));
        var z_poles = DigitalFilter.ToZArray(translated_poles, dt);

        z_poles.AssertEquals(
            (-0.19540205839647506, -0.6503947744476488),
            (-0.19540205839647506, +0.6503947744476488),
            (-0.14149393799941629, +0.19507879193594393),
            (-0.14149393799941629, -0.19507879193594393)
        );

        var Gz0 = z_poles.Multiply(z => 1 + z).Abs / (1 << N);
        Gz0.AssertEquals(0.051852987727303922);

        var B = new double[N + 1];
        for (var i = 0; i < B.Length; i++)
            B[i] = BinomialCoefficient(N, i) * (i % 2 == 0 ? Gz0 : -Gz0);
        var A = GetCoefficientsInverted(z_poles).ToRe();

        B.AssertEquals(
            +0.051852987727303922,
            -0.20741195090921569,
            +0.31111792636382352,
            -0.20741195090921569,
            +0.051852987727303922);

    }

    [TestMethod, Ignore]
    public void TransmissionCoefficient()
    {
        const double fp = 0.1;    // Гц // Граничная частота полосы пропускания
        const double fs = 0.3;    // Гц // Граничная частота полосы запирания
        const double fd = 1;      // Гц // Частота дискретизации
        const double dt = 1 / fd; // 2с // Период дискретизации
        const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
        const double Rs = 30; // Неравномерность в полосе пропускания (дБ)
        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ButterworthHighPass(dt, fp, fs, Gp, Gs);

        var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
        var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
        var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
        var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

        var transmission_0_abs = transmission_0.Abs;
        var transmission_fp_abs = transmission_fp.Abs;
        var transmission_fs_abs = transmission_fs.Abs;
        var transmission_fd05_abs = transmission_fd05.Abs;

        Assert.That.Value(transmission_0_abs).IsEqual(0, 1e-15);
        Assert.That.Value(transmission_fp_abs).LessOrEqualsThan(Gs, 1e-15);
        Assert.That.Value(transmission_fs_abs).IsEqual(Gp);
        Assert.That.Value(transmission_fd05_abs).IsEqual(1, 1e-15);
    }
}