using System;
using System.Collections;
using System.Diagnostics;
using System.Linq;
using System.Reflection;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Linq.Enumerable;
using static System.Math;
using static MathCore.SpecialFunctions.EllipticJacobi;

using static MathCore.Polynom.Array;
using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Infrastructure;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters;

[TestClassHandler("FailResultHandler")]
public class EllipticBandPass : UnitTest
{
    // ReSharper disable once UnusedMember.Local
    private static void FailResultHandler(TestResult result)
    {
        if(result.TestFailureException?.InnerException is not AssertFailedException exception) return;
        switch (exception.Data["Actual"])
        {
            case IEnumerable<Complex> actual:
                result.ToDebugEnum(actual);
                break;
            case IEnumerable actual:
                result.ToDebugEnum(actual);
                break;
            case { } actual:
                result.ToDebug(actual);
                break;
        }

        switch (exception.Data["Expected"])
        {
            case IEnumerable<Complex> expected:
                result.ToDebugEnum(expected);
                break;
            case IEnumerable expected:
                result.ToDebugEnum(expected);
                break;
            case { } expected:
                result.ToDebug(expected);
                break;
        }
    }

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

        //const double wsl = fsl * Consts.pi2; //  2
        //const double wpl = fpl * Consts.pi2; //  4
        //const double wph = fph * Consts.pi2; // 12
        //const double wsh = fsh * Consts.pi2; // 15

        //const double wc = wpl * wph; // 48
        //const double dw = wph - wpl; // 8

        //var wp = wc / wsh > wsl             // 15
        //    ? wsh
        //    : wsl;
        //var w0 = Abs(dw * wp / (wc - wp.Pow2()));
        //var f0 = w0 / Consts.pi2;
        //const double f1 = 1 / Consts.pi2;   // 0.159
        //var w1 = f1 * Consts.pi2;

        //f0.AssertEquals(0.10790165633348836);
        //w0.AssertEquals(0.67796610169491522);

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fsl = DigitalFilter.ToDigitalFrequency(fsl, dt);
        var Fpl = DigitalFilter.ToDigitalFrequency(fpl, dt);
        var Fph = DigitalFilter.ToDigitalFrequency(fph, dt);
        var Fsh = DigitalFilter.ToDigitalFrequency(fsh, dt);

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
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
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

        //var Fp = DigitalFilter.ToAnalogFrequency(f0, dt);
        //var Fs = DigitalFilter.ToAnalogFrequency(f1, dt);

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
            /*[ 0]*/ (-0.105281264621164369, 0.993710811208774136),
            /*[ 1]*/ (-0.105281264621164369, -0.993710811208774136),
            /*[ 2]*/ (-0.364290595873427325, 0.478602767640666948),
            /*[ 3]*/ (-0.364290595873427325, -0.478602767640666948)
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
        ppf_poles.AssertEquals(AccuracyComplex.Eps(1e-13),
            /*[ 0]*/ (-0.232612131820381429, -4.057810063607998785),
            /*[ 1]*/ (-0.781092257506545651, 13.625789842998447199),
            /*[ 2]*/ (-0.232612131820381429, 4.057810063607998785),
            /*[ 3]*/ (-0.781092257506545651, -13.625789842998447199),
            /*[ 4]*/ (-1.223131650509463597, -5.310820631246345513),
            /*[ 5]*/ (-2.284453268385775448, 9.919064349130303881),
            /*[ 6]*/ (-1.223131650509463597, 5.310820631246345513),
            /*[ 7]*/ (-2.284453268385775448, -9.919064349130303881)
        );

        var z_zeros_enum = DigitalFilter.ToZ(ppf_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = DigitalFilter.ToZArray(ppf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(AccuracyComplex.Eps(1e-13),
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
        z_poles.AssertEquals(AccuracyComplex.Eps(1e-13),
            /*[ 0]*/ (0.900559137768189188, -0.381172136621393376),
            /*[ 1]*/ (0.346108870590590256, 0.882619467214864506),
            /*[ 2]*/ (0.900559137768189188, 0.381172136621393376),
            /*[ 3]*/ (0.346108870590590256, -0.882619467214864506),
            /*[ 4]*/ (0.773670946459013131, -0.443838751538378096),
            /*[ 5]*/ (0.498153041879348168, 0.666845008413482043),
            /*[ 6]*/ (0.773670946459013131, 0.443838751538378096),
            /*[ 7]*/ (0.498153041879348168, -0.666845008413482043)
        );

        var Fp0 = (Fpl * Fph).Sqrt().AssertEquals(1.1853844635393842);
        var ffp0 = DigitalFilter.ToAnalogFrequency(Fp0, dt).AssertEquals(1.1347392325852204);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);
        z0.AssertEquals(new Complex(0.7564175596225313, 0.6540890424817513));

        var norm_0 = z_zeros.Multiply(z => z0 - z);
        var norm_p = z_poles.Multiply(z => z0 - z);

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * (norm_p / norm_0).Abs;
        else
            g_norm = (z0 * norm_p / norm_0).Abs;

        g_norm.AssertEquals(0.027089200894329743, 1e-16);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var h_f00 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0, dt);
        var h_fsl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsl, dt);
        var h_fpl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fpl, dt);
        var h_fcc = DoubleArrayDSPExtensions.FrequencyResponse(A, B, ffp0, dt);
        var h_fph = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fph, dt);
        var h_fsh = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsh, dt);
        var h_fd5 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fd / 2, dt);

        //var h_f00_db = h_f00.Power.In_dB_byPower().ToDebug();
        //var h_fsl_db = h_fsl.Power.In_dB_byPower().ToDebug();
        //var h_fpl_db = h_fpl.Power.In_dB_byPower().ToDebug();
        //var h_fcc_db = h_fcc.Power.In_dB_byPower().ToDebug();
        //var h_fph_db = h_fph.Power.In_dB_byPower().ToDebug();
        //var h_fsh_db = h_fsh.Power.In_dB_byPower().ToDebug();
        //var h_fd5_db = h_fd5.Power.In_dB_byPower().ToDebug();

        h_f00.Abs.AssertLessOrEqualsThan(Gs, 7e-4);
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp, 2.1e-12);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp, 4.8e-14);
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp, 0.107);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs, 4.7e-15);

        var filter = new DSP.Filters.EllipticBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);

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

    [TestMethod]
    public void TransmissionCoefficient()
    {
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

        var filter = new DSP.Filters.EllipticBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        var h_f0 = filter.FrequencyResponse(0);
        var h_sl = filter.FrequencyResponse(fsl);
        var h_pl = filter.FrequencyResponse(fpl);
        var h_c0 = filter.FrequencyResponse((fpl * fph).Sqrt());
        var h_ph = filter.FrequencyResponse(fph);
        var h_sh = filter.FrequencyResponse(fsh);
        var h_fd = filter.FrequencyResponse(fd / 2);

        //h_f0.Abs.ToDebug();
        //h_sl.Abs.ToDebug();
        //h_pl.Abs.ToDebug();
        //h_c0.Abs.ToDebug();
        //h_ph.Abs.ToDebug();
        //h_sh.Abs.ToDebug();
        //h_fd.Abs.ToDebug();

        //(-Rs).ToDebug();
        //(-Rp).ToDebug();

        var h_f0_db = h_f0.Power.In_dB_byPower();
        var h_pl_db = h_pl.Power.In_dB_byPower();
        var h_sl_db = h_sl.Power.In_dB_byPower();
        var h_c0_db = h_c0.Power.In_dB_byPower();
        var h_sh_db = h_sh.Power.In_dB_byPower();
        var h_ph_db = h_ph.Power.In_dB_byPower();
        var h_fd_db = h_fd.Power.In_dB_byPower();

        h_f0_db.AssertLessOrEqualsThan(-Rs, 2.8e-12);
        h_sl_db.AssertLessOrEqualsThan(-Rp);

        h_pl_db.AssertGreaterOrEqualsThan(-Rp, 2e-11);
        h_c0_db.AssertGreaterOrEqualsThan(-Rp);
        h_ph_db.AssertGreaterOrEqualsThan(-Rp, 3.9e-14);

        h_sh_db.AssertLessOrEqualsThan(-Rp);
        h_fd_db.AssertLessOrEqualsThan(-Rs, 4.1e-12);
    }

    [TestMethod]
    public void ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.02708967396788694, +0.05498553801823595, +0.05261996621917118, -0.02740368118281115, -0.13673170899533040,
            -0.18194040671553580, -0.09108091858457290, +0.09505079247123455, +0.23283535687203472, +0.20747423995168038,
            +0.04968683273994258, -0.09716819648256618, -0.12502844341669514, -0.05899420212080153, -0.00459152091978283,
            -0.01967283073840540, -0.06154068812320151, -0.05644509377471133, +0.00628302491495707, +0.06762607206348205,
            +0.07490722169072217, +0.03984151644980295, +0.01311970744913616, +0.01783500452400363, +0.02733953423662142,
            +0.00822095506585899, -0.03282406171369425, -0.05733700186976945, -0.04357666173802005, -0.01181631116874772,
            +0.00285976371646529, -0.00714451940809534, -0.01697171565681492, -0.00404223007325603, +0.02339915725544523,
            +0.03795344055591302, +0.02713693356164100, +0.00669318473961412, -0.00084172064188706, +0.00683463255056472,
            +0.01207965862721436, +0.00161006116246080, -0.01700967624124640, -0.02581799397591739, -0.01854345418838519,
            -0.00657465603757914, -0.00332855205344888, -0.00811514685567273, -0.00899641022394199, +0.00082786824339936,
            +0.01438907465879404, +0.01985513246192067, +0.01458675312021840, +0.00696352879189158, +0.00507225023457790,
            +0.00724533730798506, +0.00563594413757625, -0.00289583042631034, -0.01259111387523796, -0.01582100685322239,
            -0.01185346275709155, -0.00670149361643199, -0.00518270775346538, -0.00554608185316298, -0.00282467453490760,
            +0.00411680274761875, +0.01083925103554189, +0.01257897009290240, +0.00955020221600513, +0.00596439396620610,
            +0.00460845657304620, +0.00395032501698146, +0.00094533667665190, -0.00451082778558660, -0.00912542101007898,
            -0.00999318456351062, -0.00770494075006027, -0.00514510347752270, -0.00387582633024650, -0.00268893619865597,
            +0.00020310338515449, +0.00442616080236914, +0.00760107159387510, +0.00799469388704454, +0.00627547617758924,
            +0.00438079785611501, +0.00315938774350132, +0.00173200783111084, -0.00088199050434732, -0.00412973519164065,
            -0.00632558415109605, -0.00646234758087474, -0.00515333752367821, -0.00368238368381474, -0.00249971647860045,
            -0.00101026014713156, +0.00126514267201307, +0.00375183837130058, +0.00527490744468932, +0.00526872258614187
        };

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

        var filter = new DSP.Filters.EllipticBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod]
    public void SignalProcessing()
    {
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

        var filter = new DSP.Filters.EllipticBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        //const double total_time = 1 / fpl;
        //const int samples_count = (int)(total_time * fd) + 1;

        // Метод формирования гармонического сигнала на заданной частоте
        static EnumerableSignal GetSinSignal(double f0) => MathEnumerableSignal.Sin(dt, f0, (int)(100 * fd / (fsl * 0.1)));

        var x_low = GetSinSignal(fsl * 0.1);                        // частота, близкая к нулю
        var x_sl99 = GetSinSignal(fsl * 0.99);                      // частота чуть ниже нижней границы подавления
        var x_sl = GetSinSignal(fsl);                               // частота на нижней границе подавления

        var x_sl_pl = GetSinSignal(fsl + (fpl - fsl) * 0.1);        // частота выше нижней границы подавления
        var x_pl_sl = GetSinSignal(fpl - (fpl - fsl) * 0.1);        // частота ниже нижней границы пропускания

        var x_pl = GetSinSignal(fpl);                               // частота на границе пропускания
        var x_pl_ph = GetSinSignal(fpl + (fph - fpl) * 0.1);        // частота выше нижней границы пропускания
        var x_c0 = GetSinSignal((fpl * fph).Sqrt());                // частота в середине полосы пропускания
        var x_ph_pl = GetSinSignal(fph - (fph - fpl) * 0.1);        // частота ниже верхней границы пропускания
        var x_ph = GetSinSignal(fph);                               // частота на границе пропускания

        var x_ph_sh = GetSinSignal(fph + (fsh - fph) * 0.1);        // частота выше верхней границы пропускания
        var x_sh_ph = GetSinSignal(fsh - (fsh - fph) * 0.1);        // частота ниже верхней границы подавления

        var x_sh = GetSinSignal(fsh);                               // частота на верхней границе подавления
        var x_sh_fd05 = GetSinSignal(fsh + (fd / 2 - fsh) * 0.1);   // частота выше верхней границы подавления
        var x_fd05 = GetSinSignal(0.9 * (fd / 2));                  // частота ниже половины частоты дискретизации

        // Индивидуальная фильтрация каждой частотной составляющей
        var y_low = filter.ProcessIndividual(x_low);
        var y_sl99 = filter.ProcessIndividual(x_sl99);
        var y_sl = filter.ProcessIndividual(x_sl);

        var y_sl_pl = filter.ProcessIndividual(x_sl_pl);
        var y_pl_sl = filter.ProcessIndividual(x_pl_sl);

        var y_pl = filter.ProcessIndividual(x_pl);
        var y_pl_ph = filter.ProcessIndividual(x_pl_ph);
        var y_c0 = filter.ProcessIndividual(x_c0);
        var y_ph_pl = filter.ProcessIndividual(x_ph_pl);
        var y_ph = filter.ProcessIndividual(x_ph);

        var y_ph_sh = filter.ProcessIndividual(x_ph_sh);
        var y_sh_ph = filter.ProcessIndividual(x_sh_ph);

        var y_sh = filter.ProcessIndividual(x_sh);
        var y_sh_fd05 = filter.ProcessIndividual(x_sh_fd05);
        var y_fd05 = filter.ProcessIndividual(x_fd05);

        // Рассчёт коэффициентов передачи по мощности
        var k_low = y_low.Power / x_low.Power;
        var k_sl99 = y_sl99.Power / x_sl99.Power;
        var k_sl = y_sl.Power / x_sl.Power;

        var k_sl_pl = y_sl_pl.Power / x_sl_pl.Power;
        var k_pl_sl = y_pl_sl.Power / x_pl_sl.Power;

        var k_pl = y_pl.Power / x_pl.Power;
        var k_pl_ph = y_pl_ph.Power / x_pl_ph.Power;
        var k_c0 = y_c0.Power / x_c0.Power;
        var k_ph_pl = y_ph_pl.Power / x_ph_pl.Power;
        var k_ph = y_ph.Power / x_ph.Power;

        var k_ph_sh = y_ph_sh.Power / x_ph_sh.Power;
        var k_sh_ph = y_sh_ph.Power / x_sh_ph.Power;

        var k_sh = y_sh.Power / x_sh.Power;
        var k_sh_fd05 = y_sh_fd05.Power / x_sh_fd05.Power;
        var k_fd05 = y_fd05.Power / x_fd05.Power;

        // Рассчёт коэффициентов передачи по мощности в логарифмическим масштабе
        var k_low_db = k_low.In_dB_byPower();
        var k_sl99_db = k_sl99.In_dB_byPower();
        var k_sl_db = k_sl.In_dB_byPower();

        var k_sl_pl_db = k_sl_pl.In_dB_byPower();
        var k_pl_sl_db = k_pl_sl.In_dB_byPower();

        var k_pl_db = k_pl.In_dB_byPower();
        var k_pl_ph_db = k_pl_ph.In_dB_byPower();
        var k_c0_db = k_c0.In_dB_byPower();
        var k_ph_pl_db = k_ph_pl.In_dB_byPower();
        var k_ph_db = k_ph.In_dB_byPower();

        var k_ph_sh_db = k_ph_sh.In_dB_byPower();
        var k_sh_ph_db = k_sh_ph.In_dB_byPower();

        var k_sh_db = k_sh.In_dB_byPower();
        var k_sh_fd05_db = k_sh_fd05.In_dB_byPower();
        var k_fd05_db = k_fd05.In_dB_byPower();

        // Сравнение коэффициентов передачи с заданными параметрами фильтрации
        k_low_db.AssertLessOrEqualsThan(-Rs);        // Коэффициенты передачи в нижней полосе пропускания
        k_sl99_db.AssertLessOrEqualsThan(-Rs);       // должны быть меньше, чем заданный уровень
        k_sl_db.AssertLessOrEqualsThan(-Rs);         // неравномерности АЧХ (допуск) Rp - не более -1дБ

        // Коэффициенты передачи в переходной полосе
        // между полосой пропускания и полосой заграждения
        k_sl_pl_db.AssertLessOrEqualsThan(-Rp);      // Коэффициент передачи у нижнего края переходной полосы не должен быть ниже уровня подавления Rs (-40дБ), но может быть всё ещё порядка уровня пропускания Rp (-1дБ)
        k_pl_sl_db.AssertLessOrEqualsThan(-Rp);         // Коэффициент передачи у верхнего края переходной полосы должен быть гарантировано меньше коэффициента пропускания Rp (-1дБ) и должен приближаться к уровню Rs (-40дБ)

        // Коэффициенты передачи в полосе заграждения
        // должны бытьниже уровня подавления Rs
        k_pl_db.AssertGreaterOrEqualsThan(-Rs);      // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_pl_ph_db.AssertGreaterOrEqualsThan(-Rp);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertGreaterOrEqualsThan(-Rp);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_ph_pl_db.AssertGreaterOrEqualsThan(-Rp);    // Коэффициент передачи у верхнего края полосы подавления
        k_ph_db.AssertGreaterOrEqualsThan(-Rp, 3.7e-3);      // Коэффициент передачи на верхней границе полосы подавления

        // Коэффициенты передачи в переходной полосе
        // между верхней полосой заграждения и верхней
        // полосой пропускания
        k_ph_sh_db.AssertLessThan(-Rp);
        k_sh_ph_db.AssertLessOrEqualsThan(-Rp);
        k_sh_ph_db.AssertLessOrEqualsThan(-Rp);

        // Коэффициенты передачи в верхней полосе пропускания
        // от верхней границы полосы пропускания до частоты,
        // близкой к половине частоты дискретизации
        k_sh_db.AssertLessOrEqualsThan(-Rs);
        k_sh_fd05_db.AssertLessOrEqualsThan(-Rs, 0.327);
        k_fd05_db.AssertLessOrEqualsThan(-Rs);

        // Суммарный сигнал
        var x =
            x_low +
            x_sl99 +
            x_sl +
            x_sl_pl +
            x_pl_sl +
            x_pl +
            x_pl_ph +
            x_c0 +
            x_ph_pl +
            x_ph +
            x_ph_sh +
            x_sh_ph +
            x_sh +
            x_sh_fd05 +
            x_fd05;

        // Фильтрация суммарного гармонического сигнала
        var y = filter.ProcessIndividual(x);

        var X = x.GetSpectrum();    // Спектр входного сигнала фильтра
        var Y = y.GetSpectrum();    // Спектр выходного сигнала фильтра

        var H = Y / X;              // Коэффициент передачи фильтра как отношение спектров на выходе и на входе

        // Извлекаем из спетра коэффициенты передачи на частотах гармонических составляющих исходного сигнала
        var h_low = H.GetValue(fsl * 0.1).Power.In_dB_byPower();
        var h_sl99 = H.GetValue(fsl * 0.99).Power.In_dB_byPower();
        var h_sl = H.GetValue(fsl).Power.In_dB_byPower();

        var h_sl_pl = H.GetValue(fsl + (fpl - fsl) * 0.1).Power.In_dB_byPower();
        var h_pl_sl = H.GetValue(fpl - (fpl - fsl) * 0.1).Power.In_dB_byPower();

        var h_pl = H.GetValue(fpl).Power.In_dB_byPower();
        var h_pl_ph = H.GetValue(fpl + (fph - fpl) * 0.1).Power.In_dB_byPower();
        var h_c0 = H.GetValue((fpl * fph).Sqrt()).Power.In_dB_byPower();
        var h_ph_pl = H.GetValue(fph - (fph - fpl) * 0.1).Power.In_dB_byPower();
        var h_ph = H.GetValue(fph).Power.In_dB_byPower();

        var h_ph_sh = H.GetValue(fph + (fsh - fph) * 0.1).Power.In_dB_byPower();
        var h_sh_ph = H.GetValue(fsh - (fsh - fph) * 0.1).Power.In_dB_byPower();

        var h_sh = H.GetValue(fsh).Power.In_dB_byPower();
        var h_sh_fd05 = H.GetValue(fsh + (fd / 2 - fsh) * 0.1).Power.In_dB_byPower();
        var h_fd05 = H.GetValue(0.9 * (fd / 2)).Power.In_dB_byPower();

        // Тест фактических коэффициентов передачи
        h_low.AssertLessOrEqualsThan(-Rs);
        h_sl99.AssertLessOrEqualsThan(-Rs);
        h_sl.AssertLessOrEqualsThan(-Rs);

        h_sl_pl.AssertLessOrEqualsThan(-Rp);
        h_pl_sl.AssertLessOrEqualsThan(-Rp);
        h_pl_sl.AssertLessOrEqualsThan(-Rp);

        h_pl.AssertGreaterOrEqualsThan(-Rp);
        h_pl_ph.AssertGreaterOrEqualsThan(-Rp);
        h_c0.AssertGreaterOrEqualsThan(-Rp);
        h_ph_pl.AssertGreaterOrEqualsThan(-Rp);
        h_ph.AssertGreaterOrEqualsThan(-Rp, 4e-5);

        h_ph_sh.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);

        h_sh.AssertLessOrEqualsThan(-Rs);
        h_sh_fd05.AssertLessOrEqualsThan(-Rs);
        h_fd05.AssertLessOrEqualsThan(-Rs);
    }
}