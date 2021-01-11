using System;

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

            const double pi2 = 2 * Math.PI;

            const double fd = 5000;         // Частота дискретизации
            const double dt = 1 / fd;       // Период дискретизации

            const double fp = fd / pi2;     // Граничная частота конца интервала пропускания
            const double fs = 1.383 * fp;     // Граничная частота начала интервала подавления

            const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
            const double Rs = 45;   // Уровень подавления более 45 дБ

            var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
            var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления
        }
    }
}