using System;
using System.Linq;
using System.Reflection;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions.EllipticJacobi;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class EllipticBandSop : UnitTest
{
    [TestMethod]
    public void Creation()
    {
        const double fd = 10;         // Частота дискретизации
        const double dt = 1 / fd;       // Период дискретизации

        const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
        const double Rs = 40;   // Уровень подавления более 40 дБ

        var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
        var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

        const double fpl = 2 / Consts.pi2;  // нижняя частота границы полосы пропускания
        const double fsl = 4 / Consts.pi2;  // нижняя частота границы полосы заграждения
        const double fsh = 12 / Consts.pi2; // верхняя частота границы полосы заграждения
        const double fph = 15 / Consts.pi2; // верхняя частота границы полосы пропускания

        const double wpl = fpl * Consts.pi2; //  2
        const double wsl = fsl * Consts.pi2; //  4
        const double wsh = fsh * Consts.pi2; // 12
        const double wph = fph * Consts.pi2; // 15

        const double wc = wsl * wsh; // 45
        const double dw = wsh - wsl; // 4

        var wp = wc / wph > wpl             // 3
            ? wph
            : wpl;
        var w0 = Abs(dw * wp / (wc - wp.Pow2()));
        var f0 = w0 / Consts.pi2;
        const double f1 = 1 / Consts.pi2;   // 0.159
        var w1 = f1 * Consts.pi2;

        f0.AssertEquals(0.10790165633348836);
        w0.AssertEquals(0.67796610169491522);

        var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.50884713990958752);
        eps_s.AssertEquals(99.994999874993752);

        var kw = f0 / f1;
        var k_eps = eps_p / eps_s;

        var Kw = FullEllipticIntegral(kw);
        var Tw = FullEllipticIntegralComplimentary(kw);
        var K_eps = FullEllipticIntegral(k_eps);
        var T_eps = FullEllipticIntegralComplimentary(k_eps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N);

        double_N.AssertEquals(4.0903829785068888);
        N.AssertEquals(5);

        var Fpl = DigitalFilter.ToAnalogFrequency(fpl, dt);
        var Fsl = DigitalFilter.ToAnalogFrequency(fsl, dt);
        var Fsh = DigitalFilter.ToAnalogFrequency(fsh, dt);
        var Fph = DigitalFilter.ToAnalogFrequency(fph, dt);

        Fpl.AssertEquals(0.31830989042792263);
        Fsl.AssertEquals(0.63661980632063808);
        Fsh.AssertEquals(1.9098602338357447);
        Fph.AssertEquals(2.3873259368731512);

        var L = N / 2;
        var r = N % 2;

        var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);
        var m = (1 - k_eps * k_eps).Sqrt();

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));

        var k_W = Sqrt(1 - kp * kp);

        var im_pz = Enumerable.Range(0, L).ToArray(i => 1 / (k_W * cd_uk(u[i], k_W)));

        var v0_complex = sn_inverse((0, 1 / eps_p), k_eps) / N;

        var zeros = new Complex[N - r];
        var poles = new Complex[N];

        if (r != 0) poles[0] = Complex.ImValue(w0) * sn_uk(v0_complex, k_W);
        for (var i = 0; i < L; i++)
        {
            // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_W);

            var pp = new Complex(-p_re * w0, p_im * w0);

            poles[r + 2 * i] = pp;
            poles[r + 2 * i + 1] = pp.ComplexConjugate;

            var p0_im = w0 / (k_W * cd_uk(u[i], k_W));
            zeros[2 * i] = (0, p0_im);
            zeros[2 * i + 1] = zeros[2 * i].ComplexConjugate;
        }

        zeros.AssertEquals(
            (0, +0.85003902981662161),
            (0, -0.85003902981662161),
            (0, +1.1961277565481045),
            (0, -1.1961277565481045)
            );

        poles.AssertEquals(
            -0.26125040018686613,
            (-0.033844548398070828, +0.67674444106992981),
            (-0.033844548398070828, -0.67674444106992981),
            (-0.14854693514997036, +0.5023959058648525),
            (-0.14854693514997036, -0.5023959058648525)
            );

        var pzf_zeros = AnalogBasedFilter.TransformToBandStop(zeros, Fsl, Fsh);
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh);

        var wc_sqrt = wc.Sqrt();
        if (zeros.Length != poles.Length)
            pzf_zeros = pzf_zeros.AppendLast(
                (0, +wc_sqrt),
                (0, -wc_sqrt));

        pzf_zeros.AssertEquals(
            (+5.128308140332369E-16, +3.6694932282406674),
            (-5.128308140332369E-16, -13.08083231510096),
            (+5.128308140332369E-16, +13.08083231510096),
            (-5.128308140332369E-16, -3.6694932282406674),

            (+4.7106421190832123E-16, +4.3489356903400047),
            (-4.7106421190832123E-16, -11.037189100458273),
            (+4.7106421190832123E-16, +11.037189100458273),
            (-4.7106421190832123E-16, -4.3489356903400047),
            (0, +wc_sqrt),
            (0, -wc_sqrt)
            );

        pzf_poles.AssertEquals(
            -1.6571848295652529,
            -28.964799063847966,

            (-0.10370546264251768, +3.1986700805557033),
            (-0.48601305764408786, -14.99048735363515),
            (-0.10370546264251768, -3.1986700805557033),
            (-0.48601305764408786, +14.99048735363515),

            (-0.57541229376599867, +2.6505995694839637),
            (-3.75433472888858, -17.294100463097848),
            (-0.57541229376599867, -2.6505995694839637),
            (-3.75433472888858, +17.294100463097848)
        );

        var h_f00 = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            0)
           .Abs;

        var h_hpl = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            fpl)
           .Abs;

        var h_hsl = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            fsl)
            .Abs;

        var fc = (fsl * fsh).Sqrt();
        var h_hfc = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            fc)
           .Abs;

        var h_hsh = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            fsh).Abs;

        var h_hph = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            fph)
           .Abs;

        //var filter = new DSP.Filters.EllipticBandStop(dt, fsl, fpl, fph, fsh, Gp, Gs);

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

        Assert.Fail();
    }
}