using System;

using MathCore.DSP.Tests.Service;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class EllipticBandPass : UnitTest
    {
        [TestMethod]
        public void Creation()
        {
            // http://www.dsplib.ru/content/filters/ch8/ch8.html

            const double fd = 100;         // Частота дискретизации
            const double dt = 1 / fd;       // Период дискретизации

            const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
            const double Rs = 40;   // Уровень подавления более 45 дБ

            var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
            var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

            const double fsl = 4;
            const double fpl = 5;
            const double fph = 15;
            const double fsh = 18;

            var filter = new DSP.Filters.EllipticBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

            var h_f0 = filter.GetTransmissionCoefficient(0).Power.In_dB_byPower();
            var h_sl = filter.GetTransmissionCoefficient(fsl).Power.In_dB_byPower();
            var h_pl = filter.GetTransmissionCoefficient(fpl).Power.In_dB_byPower();
            var h_c0 = filter.GetTransmissionCoefficient((fpl * fph).Sqrt()).Power.In_dB_byPower();
            var h_ph = filter.GetTransmissionCoefficient(fph).Power.In_dB_byPower();
            var h_sh = filter.GetTransmissionCoefficient(fsh).Power.In_dB_byPower();
            var h_fd = filter.GetTransmissionCoefficient(fd/2).Power.In_dB_byPower();

            h_f0.AssertThanValue().LessThan(-Rs);
            h_sl.AssertThanValue().LessThan(-Rs);

            h_pl.AssertThanValue().GreaterThan(-Rp);
            h_c0.AssertThanValue().GreaterThan(-Rp);
            h_ph.AssertThanValue().GreaterThan(-Rp);

            h_sh.AssertThanValue().LessThan(-Rs);
            h_fd.AssertThanValue().LessThan(-Rs);
        }
    }
}