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

        const double wc = wsl * wsh; // 48
        const double dw = wsh - wsl; // 8

        var wp = wc / wph > wpl             // 15
            ? wph
            : wpl;
        var w0 = Abs(dw * wp / (wc - wp.Pow2()));
        var f0 = w0 / Consts.pi2;
        const double f1 = 1 / Consts.pi2;   // 0.159
        var w1 = f1 * Consts.pi2;

        f0.AssertEquals(0.10790165633348836);
        w0.AssertEquals(0.67796610169491522);

        var Fpl = DigitalFilter.ToAnalogFrequency(fpl, dt);
        var Fsl = DigitalFilter.ToAnalogFrequency(fsl, dt);
        var Fsh = DigitalFilter.ToAnalogFrequency(fsh, dt);
        var Fph = DigitalFilter.ToAnalogFrequency(fph, dt);

        Fpl.AssertEquals(0.31937518051807723);
        Fsl.AssertEquals(0.64524608331077715);
        Fsh.AssertEquals(2.1776750959738593);
        Fph.AssertEquals(2.9653636313402);

        var Wpl = Consts.pi2 * Fpl;
        var Wsl = Consts.pi2 * Fsl;
        var Wsh = Consts.pi2 * Fsh;
        var Wph = Consts.pi2 * Fph;

        var Wc = Wsl * Wsh;
        var dW = Wsh - Wsl;

        var Wp = Wc / Wph > Wpl
            ? Wph
            : Wpl;
        var W0 = Abs(dW * Wp / (Wc - Wp.Pow2()));
        const double W1 = 1;
        var F0 = W0 / Consts.pi2;
        const double F1 = 1 / Consts.pi2;

        var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.50884713990958752);
        eps_s.AssertEquals(99.994999874993752);

        var kw = W0 / W1;
        var k_eps = eps_p / eps_s;

        kw.AssertEquals(0.615059351152204);

        var Kw = FullEllipticIntegral(kw);
        var Tw = FullEllipticIntegralComplimentary(kw);
        var K_eps = FullEllipticIntegral(k_eps);
        var T_eps = FullEllipticIntegralComplimentary(k_eps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N);

        double_N.AssertEquals(3.7906792606389264);
        N.AssertEquals(4);

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

        if (r != 0) poles[0] = Complex.ImValue(W0) * sn_uk(v0_complex, k_W);
        for (var i = 0; i < L; i++)
        {
            // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_W);

            var pp = new Complex(-p_re * W0, p_im * W0);

            poles[r + 2 * i] = pp;
            poles[r + 2 * i + 1] = pp.ComplexConjugate;

            var p0_im = W0 / (k_W * cd_uk(u[i], k_W));
            zeros[2 * i] = (0, p0_im);
            zeros[2 * i + 1] = zeros[2 * i].ComplexConjugate;
        }

        zeros.AssertEquals(
            (0, +0.98996902542422383),
            (0, -0.98996902542422383),
            (0, +2.1682610011632977),
            (0, -2.1682610011632977)
            );

        poles.AssertEquals(
            (-0.064754226306376769, +0.61119112677499909),
            (-0.064754226306376769, -0.61119112677499909),
            (-0.22406033752875973, +0.29436910772471786),
            (-0.22406033752875973, -0.29436910772471786)
            );

        var pzf_zeros = AnalogBasedFilter.TransformToBandStop(zeros, Fsl, Fsh);
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh);

        var wc_sqrt = wc.Sqrt();
        if (zeros.Length != poles.Length)
            pzf_zeros = pzf_zeros.AppendLast(
                (0, +wc_sqrt),
                (0, -wc_sqrt));

        pzf_zeros.AssertEquals(
            (+5.44664340994985E-16, +4.0319948722574663),
            (-5.44664340994985E-16, -13.758092567621574),
            (+5.44664340994985E-16, +13.758092567621574),
            (-5.44664340994985E-16, -4.0319948722574663),

            (+4.758917038943161E-16, +5.551565428550191),
            (-4.758917038943161E-16, -9.9922372164459166),
            (+4.758917038943161E-16, +9.9922372164459166),
            (-4.758917038943161E-16, -5.551565428550191)
            );

        pzf_poles.AssertEquals(
            (-0.22795541315219781, +2.9727034750401131),
            (-1.4225863599870872, -18.551555137033667),
            (-0.22795541315219781, -2.9727034750401131),
            (-1.4225863599870872, +18.551555137033667),

            (-1.1307190536141292, +1.7343356882717966),
            (-14.633073912059702, -22.444710941843681),
            (-1.1307190536141292, -1.7343356882717966),
            (-14.633073912059702, +22.444710941843681)
        );

        var h_F00 = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            0)
           .Abs;

        var h_Fpl = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fpl)
           .Abs;

        var h_Fsl = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fsl)
            .Abs;

        var Fc = (Fsl * Fsh).Sqrt();
        var h_Fc = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fc)
           .Abs;

        var h_Fsh = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fsh).Abs;

        var h_Fph = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fph)
           .Abs;

        h_F00.AssertGreaterOrEqualsThan(Gp);
        h_Fpl.AssertGreaterOrEqualsThan(Gp);
        h_Fsl.AssertLessOrEqualsThan(Gs);
        h_Fc.AssertLessOrEqualsThan(Gs, 1e-2);
        h_Fsh.AssertLessOrEqualsThan(Gs, 1e-3);
        h_Fph.AssertGreaterOrEqualsThan(Gp);

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        z_zeros.AssertEquals(
            (0.92188968196320531, +0.38745246713600884),
            (0.35757714717701156, -0.9338836029274471),
            (0.35757714717701156, +0.9338836029274471),
            (0.92188968196320531, -0.38745246713600884),

            (0.85692452818129894, +0.51544190070390872),
            (0.60049677950867342, -0.79962717425042007),
            (0.60049677950867342, +0.79962717425042007),
            (0.85692452818129894, -0.51544190070390872));

        z_poles.AssertEquals(
            (0.93565642115019709, +0.28446436884548237),
            (0.067011448261103043, -0.92401175944069935),
            (0.93565642115019709, -0.28446436884548237),
            (0.067011448261103043, +0.92401175944069935),

            (+0.88031182726988344, +0.15432943378971034),
            (-0.18664227822550081, -0.52711402412329078),
            (+0.88031182726988344, -0.15432943378971034),
            (-0.18664227822550081, +0.52711402412329078)
            );

        var G_norm = (N.IsOdd() ? 1 : Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        Assert.Fail();

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
    }
}