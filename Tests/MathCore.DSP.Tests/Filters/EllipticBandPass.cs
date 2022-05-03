using System;
using System.Diagnostics;
using System.Linq;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Linq.Enumerable;
using static System.Math;
using static MathCore.SpecialFunctions.EllipticJacobi;

using static MathCore.Polynom.Array;
using MathCore.DSP.Signals;
using static MathCore.SpecialFunctions.Distribution;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class EllipticBandPass : UnitTest
{
    [TestMethod]
    public void Creation()
    {
        // http://www.dsplib.ru/content/filters/ch8/ch8.html

        const double fd = 10;           // Частота дискретизации
        const double dt = 1 / fd;       // Период дискретизации

        const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
        const double Rs = 40;   // Уровень подавления более 40 дБ

        var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
        var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

        const double fsl = 2 / Consts.pi2;  // нижняя частота границы полосы пропускания
        const double fpl = 4 / Consts.pi2;  // нижняя частота границы полосы заграждения
        const double fph = 12 / Consts.pi2; // верхняя частота границы полосы заграждения
        const double fsh = 15 / Consts.pi2; // верхняя частота границы полосы пропускания

        const double wsl = fsl * Consts.pi2; //  2
        const double wpl = fpl * Consts.pi2; //  4
        const double wph = fph * Consts.pi2; // 12
        const double wsh = fsh * Consts.pi2; // 15

        const double wc = wpl * wph; // 48
        const double dw = wph - wpl; // 8

        var wp = wc / wsh > wsl             // 15
            ? wsh
            : wsl;
        var w0 = Abs(dw * wp / (wc - wp.Pow2()));
        var f0 = w0 / Consts.pi2;
        const double f1 = 1 / Consts.pi2;   // 0.159
        var w1 = f1 * Consts.pi2;

        f0.AssertEquals(0.10790165633348836);
        w0.AssertEquals(0.67796610169491522);

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fsl = DigitalFilter.ToAnalogFrequency(fsl, dt);
        var Fpl = DigitalFilter.ToAnalogFrequency(fpl, dt);
        var Fph = DigitalFilter.ToAnalogFrequency(fph, dt);
        var Fsh = DigitalFilter.ToAnalogFrequency(fsh, dt);

        Fsl.AssertEquals(0.31937518051807723);
        Fpl.AssertEquals(0.64524608331077715);
        Fph.AssertEquals(2.1776750959738593);
        Fsh.AssertEquals(2.9653636313402);

        var Wsl = Consts.pi2 * Fsl;
        var Wpl = Consts.pi2 * Fpl;
        var Wph = Consts.pi2 * Fph;
        var Wsh = Consts.pi2 * Fsh;

        //Wsl.ToDebug();
        //Wpl.ToDebug();
        //Wph.ToDebug();
        //Wsh.ToDebug();

        var Wc = Wpl * Wph;
        var dW = Wph - Wpl;

        // Выбор опорной частоты
        // Если   Wc / Wsh > Wsl
        // то есть      Wc > Wsl*Wsh
        // то есть Wpl*Wph > Wsl*Wsh
        // то есть центральная частота по границам подавления > центральной частоты по границам пропускания
        // то выбираем в качестве опорной частоты выбираем верхнюю границу пропускания
        // иначе, выбираем нижнюю границу пропускания
        var Wp = Wc / Wsh > Wsl
            ? Wsh
            : Wsl;
        //Wp.ToDebug();
        var W0 = Abs(dW * Wp / (Wc - Wp.Pow2()));   // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F0 = W0 / Consts.pi2;
        const double F1 = 1 / Consts.pi2;

        W0.AssertEquals(0.615059351152204);
        F0.AssertEquals(0.0978897360307671);

        var eps_p = (1 / (Gp * Gp) - 1).Sqrt();
        //var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.5088471399095873);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_p / eps_s;
        var kW = F0 / F1;

        kEps.AssertEquals(0.005088725841749186);
        kW.AssertEquals(0.615059351152204);

        var Kw = FullEllipticIntegral(kW);
        var Tw = FullEllipticIntegralComplimentary(kW);
        var K_eps = FullEllipticIntegral(kEps);
        var T_eps = FullEllipticIntegralComplimentary(kEps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N);
        N.AssertEquals(4);

        var Fp = DigitalFilter.ToAnalogFrequency(f0, dt);
        var Fs = DigitalFilter.ToAnalogFrequency(f1, dt);

        var (L, r) = N.GetDivMod(2);
        (L, r).AssertEquals((2, 0));

        var u = Range(1, L).ToArray(i => (2 * i - 1d) / N);
        var m = (1 - kEps * kEps).Sqrt();

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));

        var k_W = Sqrt(1 - kp * kp);

        //var im_pz = Range(0, L).ToArray(i => 1 / (k_W * cd_uk(u[i], k_W)));

        var v0_complex = sn_inverse((0, 1 / eps_p), kEps) / N;

        var zeros = new Complex[N - r];
        var poles = new Complex[N];

        if (r != 0) poles[0] = Complex.i * sn_uk(v0_complex, k_W);
        for (var i = 0; i < L; i++)
        {
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_W);
            (poles[r + 2 * i], poles[r + 2 * i + 1]) = Complex.Conjugate(-p_re, p_im);

            var p0_im = 1 / (k_W * cd_uk(u[i], k_W));
            (zeros[2 * i], zeros[2 * i + 1]) = Complex.Conjugate(0, p0_im);
        }

        zeros.ToRe().Sum(v => v * v).AssertEquals(0);
        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ +1.6095504012250093,
            /*[ 1]*/ -1.6095504012250093,
            /*[ 2]*/ +3.5252874329956083,
            /*[ 3]*/ -3.5252874329956083
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.105281264621164411, 0.993710811208774136),
            /*[ 1]*/ (-0.105281264621164411, -0.993710811208774136),
            /*[ 2]*/ (-0.364290595873427714, 0.478602767640667393),
            /*[ 3]*/ (-0.364290595873427714, -0.478602767640667393)
        );

        var ppf_zeros = AnalogBasedFilter.TransformToBandPassW(zeros, Wpl, Wph).ToArray();
        var ppf_poles = AnalogBasedFilter.TransformToBandPassW(poles, Wpl, Wph).ToArray();

        ppf_zeros.ToRe().Sum(v => v * v).AssertEquals(0, 7e-30);
        //ppf_zeros.ToIm().ToDebugEnum();
        ppf_zeros.ToIm().AssertEquals(
            /*[ 0]*/ +18.4966696758644030,
            /*[ 1]*/ -02.9990565683874317,
            /*[ 2]*/ +02.9990565683874317,
            /*[ 3]*/ -18.4966696758644030,
            /*[ 4]*/ +35.5057109896055700,
            /*[ 5]*/ -01.5623559460880400,
            /*[ 6]*/ +01.5623559460880400,
            /*[ 7]*/ -35.5057109896055700
        );

        //ppf_poles.ToDebugEnum();
        ppf_poles.AssertEquals(
            /*[ 0]*/ (-0.232612131820379653, -04.057810063607998785),
            /*[ 1]*/ (-0.781092257506547871, +13.625789842998447199),
            /*[ 2]*/ (-0.232612131820379653, +04.057810063607998785),
            /*[ 3]*/ (-0.781092257506547871, -13.625789842998447199),
            /*[ 4]*/ (-1.223131650509463597, -5.3108206312463428490),
            /*[ 5]*/ (-2.284453268385779001, +09.919064349130305658),
            /*[ 6]*/ (-1.223131650509463597, +05.310820631246342849),
            /*[ 7]*/ (-2.284453268385779001, -09.919064349130305658)
        );

        var z_zeros_enum = DigitalFilter.ToZ(ppf_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = DigitalFilter.ToZArray(ppf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 1]*/ (+0.077982915792994864, +0.996954695482409003),
            /*[ 2]*/ (+0.956017287213403399, -0.293310324654836141),
            /*[ 3]*/ (+0.956017287213403399, +0.293310324654836141),
            /*[ 4]*/ (+0.077982915792994864, -0.996954695482409003),
            /*[ 5]*/ (-0.518262521157172529, +0.855221584832732806),
            /*[ 6]*/ (+0.987869246083113217, -0.155287966833175056),
            /*[ 7]*/ (+0.987869246083113217, +0.155287966833175056),
            /*[ 8]*/ (-0.518262521157172529, -0.855221584832732806)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.900559137768189522, -0.381172136621393431),
            /*[ 1]*/ (0.346108870590590256, 0.882619467214864506),
            /*[ 2]*/ (0.900559137768189522, 0.381172136621393431),
            /*[ 3]*/ (0.346108870590590256, -0.882619467214864506),
            /*[ 4]*/ (0.773670946459013353, -0.443838751538377929),
            /*[ 5]*/ (0.498153041879348002, 0.666845008413482154),
            /*[ 6]*/ (0.773670946459013353, 0.443838751538377929),
            /*[ 7]*/ (0.498153041879348002, -0.666845008413482154)
        );

        var Fp0 = (Fpl * Fph).Sqrt().AssertEquals(1.1853844635393842);
        var ffp0 = DigitalFilter.ToDigitalFrequency(Fp0, dt).AssertEquals(1.1347392325852204);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);
        z0.AssertEquals(new Complex(0.7564175596225313, 0.6540890424817513));

        var norm_0 = z_zeros.Multiply(z => z0 - z);
        var norm_p = z_poles.Multiply(z => z0 - z);

        double G_norm;
        if (N.IsEven())
            G_norm = Gp * (norm_p / norm_0).Abs;
        else
            G_norm = (z0 * norm_p / norm_0).Abs;

        G_norm.AssertEquals(0.027089200894329788);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var h_f00 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0, dt);
        var h_fsl = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fsl, dt);
        var h_fpl = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fpl, dt);
        var h_fcc = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, ffp0, dt);
        var h_fph = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fph, dt);
        var h_fsh = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fsh, dt);
        var h_fd5 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fd / 2, dt);

        //var h_f00_db = h_f00.Power.In_dB_byPower().ToDebug();
        //var h_fsl_db = h_fsl.Power.In_dB_byPower().ToDebug();
        //var h_fpl_db = h_fpl.Power.In_dB_byPower().ToDebug();
        //var h_fcc_db = h_fcc.Power.In_dB_byPower().ToDebug();
        //var h_fph_db = h_fph.Power.In_dB_byPower().ToDebug();
        //var h_fsh_db = h_fsh.Power.In_dB_byPower().ToDebug();
        //var h_fd5_db = h_fd5.Power.In_dB_byPower().ToDebug();

        h_f00.Abs.AssertLessOrEqualsThan(Gs, 7e-4);
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp, 3.1e-14);
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp, 0.107);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs, 4.7e-15);

        //var filter = new DSP.Filters.EllipticBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        //var h_f0 = filter.GetTransmissionCoefficient(0).Power.In_dB_byPower();
        //var h_sl = filter.GetTransmissionCoefficient(fsl).Power.In_dB_byPower();
        //var h_pl = filter.GetTransmissionCoefficient(fpl).Power.In_dB_byPower();
        //var h_c0 = filter.GetTransmissionCoefficient((fpl * fph).Sqrt()).Power.In_dB_byPower();
        //var h_ph = filter.GetTransmissionCoefficient(fph).Power.In_dB_byPower();
        //var h_sh = filter.GetTransmissionCoefficient(fsh).Power.In_dB_byPower();
        //var h_fd = filter.GetTransmissionCoefficient(fd/2).Power.In_dB_byPower();

        //h_f0.AssertThatValue().LessThan(-Rs);
        //h_sl.AssertThatValue().LessThan(-Rs);

        //h_pl.AssertThatValue().GreaterThan(-Rp, 1.1);
        //h_c0.AssertThatValue().GreaterThan(-Rp);
        //h_ph.AssertThatValue().GreaterThan(-Rp);

        //h_sh.AssertThatValue().LessThan(-Rs);
        //h_fd.AssertThatValue().LessThan(-Rs);
    }
}