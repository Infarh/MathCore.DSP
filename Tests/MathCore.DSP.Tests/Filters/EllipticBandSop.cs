using System;
using System.Linq;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class EllipticBandSop : UnitTest
{
    [TestMethod]
    public void Creation()
    {
        const double fd = 5000;         // Частота дискретизации
        const double dt = 1 / fd;       // Период дискретизации

        const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
        const double Rs = 40;   // Уровень подавления более 40 дБ

        var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
        var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

        const double fpl = 3 / Consts.pi2;
        const double fsl = 5 / Consts.pi2;
        const double fsh = 15 / Consts.pi2;
        const double fph = 20 / Consts.pi2;

        const double wpl = fpl * Consts.pi2; //  3
        const double wsl = fsl * Consts.pi2; //  5
        const double wsh = fsh * Consts.pi2; // 15
        const double wph = fph * Consts.pi2; // 20

        const double wc = wpl * wph; // 60
        const double dw = wph - wpl; // 17

        const double f0 = 1 / Consts.pi2;
        var ws = wc / wph > wpl
            ? wph
            : wpl;
        var f1 = Math.Abs((wc - ws.Pow2()) / (dw * ws)) / Consts.pi2;

        var eps_p = Math.Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Math.Sqrt(10.Power(Rs / 10) - 1);

        var kw = f1 / f0;
        var k_eps = eps_p / eps_s;

        var Kw = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(kw);
        var Tw = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(kw);
        var K_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(k_eps);
        var T_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(k_eps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Math.Ceiling(double_N);

        var Fp = DigitalFilter.ToAnalogFrequency(f0, dt);
        var Fs = DigitalFilter.ToAnalogFrequency(f1, dt);

        var L = N / 2;
        var r = N % 2;

        var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);
        var m = (1 - k_eps * k_eps).Sqrt();

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * SpecialFunctions.EllipticJacobi.sn_uk(ui, m).Power(4));

        var k_W = Math.Sqrt(1 - kp * kp);

        var im_pz = Enumerable.Range(0, L).ToArray(i => 1 / (k_W * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k_W)));

        var v0_complex = SpecialFunctions.EllipticJacobi.sn_inverse((0, 1 / eps_p), k_eps) / N;

        var zeros = new Complex[N - r];
        var poles = new Complex[N];

        //if (r != 0) poles[0] = Complex.i * SpecialFunctions.EllipticJacobi.sn_uk(v0_complex, k_W);
        //for (var i = 0; i < L; i++)
        //{
        //    var (p_im, p_re) = SpecialFunctions.EllipticJacobi.cd_uk(u[i] - v0_complex, k_W);

        //    poles[r + 2 * i] = (-p_re, p_im);
        //    poles[r + 2 * i + 1] = poles[r + 2 * i].ComplexConjugate;

        //    var p0_im = 1 / (k_W * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k_W));
        //    zeros[2 * i] = (0, p0_im);
        //    zeros[2 * i + 1] = zeros[2 * i].ComplexConjugate;
        //}


        //var filter = new DSP.Filters.EllipticBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        //var h_f0 = filter.GetTransmissionCoefficient(0).Power.In_dB_byPower();
        //var h_sl = filter.GetTransmissionCoefficient(fsl).Power.In_dB_byPower();
        //var h_pl = filter.GetTransmissionCoefficient(fpl).Power.In_dB_byPower();
        //var h_c0 = filter.GetTransmissionCoefficient((fpl * fph).Sqrt()).Power.In_dB_byPower();
        //var h_ph = filter.GetTransmissionCoefficient(fph).Power.In_dB_byPower();
        //var h_sh = filter.GetTransmissionCoefficient(fsh).Power.In_dB_byPower();
        //var h_fd = filter.GetTransmissionCoefficient(fd / 2).Power.In_dB_byPower();

        //h_f0.AssertThatValue().LessThan(-Rs);
        //h_sl.AssertThatValue().LessThan(-Rs);

        //h_pl.AssertThatValue().GreaterThan(-Rp, 1.1);
        //h_c0.AssertThatValue().GreaterThan(-Rp);
        //h_ph.AssertThatValue().GreaterThan(-Rp);

        //h_sh.AssertThatValue().LessThan(-Rs);
        //h_fd.AssertThatValue().LessThan(-Rs);
    }
}