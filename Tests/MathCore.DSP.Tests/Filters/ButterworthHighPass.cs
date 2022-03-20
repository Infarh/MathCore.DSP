using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ButterworthHighPass
{
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