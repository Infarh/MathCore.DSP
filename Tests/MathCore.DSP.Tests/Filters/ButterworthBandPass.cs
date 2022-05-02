using MathCore.DSP.Filters;
using System;
using System.Diagnostics;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;

using static MathCore.Polynom.Array;
using MathCore.DSP.Signals;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ButterworthBandPass
{
    [TestMethod]
    public void Creation()
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

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.5128464567198);
        kW.AssertEquals(1.6258593550796006);

        var N = (int)Ceiling(Log(kEps) / Log(kW));
        N.AssertEquals(11); // порядок фильтра - число полюсов

        //var L = N / 2;  // число парных полюсов
        var r = N % 2;  // число (1 или 0) непарных (чисто действительных) полюсов

        var alpha = eps_p.Pow(-1d / N);
        alpha.AssertEquals(1.0633442289443011);

        var poles = new Complex[N];
        if (r != 0) poles[0] = -alpha;
        for (var (i, th0) = (r, Consts.pi05 / N); i < N; i += 2)
        {
            var w = th0 * (i - r + 1);
            var sin = -alpha * Sin(w);
            var cos = alpha * Cos(w);
            (poles[i], poles[i + 1]) = Complex.Conjugate(sin, cos);
        }

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/  -1.063344228944301140,
            /*[ 1]*/ (-0.151329661971039309, 1.052520917909416820),
            /*[ 2]*/ (-0.151329661971039309, -1.052520917909416820),
            /*[ 3]*/ (-0.441729156692377722, 0.967251932723316554),
            /*[ 4]*/ (-0.441729156692377722, -0.967251932723316554),
            /*[ 5]*/ (-0.696342382202948307, 0.803621948416712528),
            /*[ 6]*/ (-0.696342382202948307, -0.803621948416712528),
            /*[ 7]*/ (-0.894542089215041636, 0.574887293173139136),
            /*[ 8]*/ (-0.894542089215041636, -0.574887293173139136),
            /*[ 9]*/ (-1.020271316205582313, 0.299578688423056627),
            /*[10]*/ (-1.020271316205582313, -0.299578688423056627)
        );

        //(Fpl, Fph).ToDebug();
        var ppf_poles = AnalogBasedFilter.TransformToBandPassW(poles, Wpl, Wph).ToArray();

        //pzf_poles.ToDebugEnum();
        ppf_poles.AssertEquals(
            /*[ 0]*/ (-5.119223805512707948, +05.409815737505828892),
            /*[ 1]*/ (-5.119223805512707948, -05.409815737505828892),
            /*[ 2]*/ (-0.317820611548505549, -03.920994507520513928),
            /*[ 3]*/ (-1.139262404384079108, +14.055229484488080516),
            /*[ 4]*/ (-0.317820611548505549, +03.920994507520513928),
            /*[ 5]*/ (-1.139262404384079108, -14.055229484488080516),
            /*[ 6]*/ (-0.975119358508945933, -03.943393058840062082),
            /*[ 7]*/ (-3.278085488944309400, +13.256612588589813839),
            /*[ 8]*/ (-0.975119358508945933, +03.943393058840062082),
            /*[ 9]*/ (-3.278085488944309400, -13.256612588589813839),
            /*[10]*/ (-1.704126634532570517, -04.000002608568959950),
            /*[11]*/ (-5.000630682483883760, +11.737705032649788706),
            /*[12]*/ (-1.704126634532570517, +04.000002608568959950),
            /*[13]*/ (-5.000630682483883760, -11.737705032649788706),
            /*[14]*/ (-2.578473111862299216, -04.129606731608852499),
            /*[15]*/ (-6.034657111619797121, +09.664929417509943832),
            /*[16]*/ (-2.578473111862299216, +04.129606731608852499),
            /*[17]*/ (-6.034657111619797121, -09.664929417509943832),
            /*[18]*/ (-3.712116727017849449, -04.462463918595362067),
            /*[19]*/ (-6.111601816481171845, +07.346967942136580554),
            /*[20]*/ (-3.712116727017849449, +04.462463918595362067),
            /*[21]*/ (-6.111601816481171845, -07.346967942136580554)
        );


        var z_zeros = Enumerable.Range(0, N)
           .Select(_ => new Complex(1))
           .Concat(Enumerable.Range(0, N).Select(_ => new Complex(-1)))
           .ToArray();
        var z_poles = DigitalFilter.ToZArray(ppf_poles, dt);

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.521820327296527386, 0.327747689180517110),
            /*[ 1]*/ (0.521820327296527386, -0.327747689180517110),
            /*[ 2]*/ (0.898027882838065006, -0.366287164652825625),
            /*[ 3]*/ (0.312146214691436663, 0.872429501649296224),
            /*[ 4]*/ (0.898027882838065006, 0.366287164652825625),
            /*[ 5]*/ (0.312146214691436663, -0.872429501649296224),
            /*[ 6]*/ (0.841918239241153277, -0.346286830383539024),
            /*[ 7]*/ (0.297539413493185079, 0.738934365168244955),
            /*[ 8]*/ (0.841918239241153277, 0.346286830383539024),
            /*[ 9]*/ (0.297539413493185079, -0.738934365168244955),
            /*[10]*/ (0.782426774908169920, -0.328495675926993413),
            /*[11]*/ (0.310983664721538322, 0.615502054910353480),
            /*[12]*/ (0.782426774908169920, 0.328495675926993413),
            /*[13]*/ (0.310983664721538322, -0.615502054910353480),
            /*[14]*/ (0.714252978329938615, -0.313537173391612911),
            /*[15]*/ (0.350320367340023753, 0.501283768993547207),
            /*[16]*/ (0.714252978329938615, 0.313537173391612911),
            /*[17]*/ (0.350320367340023753, -0.501283768993547207),
            /*[18]*/ (0.629200287496307298, -0.306604744857453426),
            /*[19]*/ (0.419506704928794050, 0.399403695263809477),
            /*[20]*/ (0.629200287496307298, 0.306604744857453426),
            /*[21]*/ (0.419506704928794050, -0.399403695263809477)
        );

        var g_norm = Complex.Real;
        for (var i = 0; i < N; i++)
        {
            var (p1, p2) = (ppf_poles[2 * i], ppf_poles[2 * i + 1]);
            var p = (2 * fd - p1) * (2 * fd - p2);
            g_norm *= 2 * fd * dW / p;
        }

        var G_norm = g_norm / eps_p;

        G_norm.Re.AssertEquals(6.892053043535363E-06);
        G_norm.Im.Abs().AssertLessThan(1e-20);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var h_f00 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0, dt);
        var h_fsl = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fsl, dt);
        var h_fpl = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fpl, dt);
        var h_fcc = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, (fsl * fsh).Sqrt(), dt);
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

        h_f00.Abs.AssertLessOrEqualsThan(Gs);
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp, 1.5e-6);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs);

        var filter = new DSP.Filters.ButterworthBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);

        //filter.B.Zip(B, (fb, b) => (fb - b).Abs()).Sum().ToDebug();
        //filter.A.Zip(A, (fa, a) => (fa - a).Abs()).Sum().ToDebug();
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

        var filter = new DSP.Filters.ButterworthBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        var h_f0 = filter.GetTransmissionCoefficient(0);
        var h_sl = filter.GetTransmissionCoefficient(fsl);
        var h_pl = filter.GetTransmissionCoefficient(fpl);
        var h_c0 = filter.GetTransmissionCoefficient((fpl * fph).Sqrt());
        var h_ph = filter.GetTransmissionCoefficient(fph);
        var h_sh = filter.GetTransmissionCoefficient(fsh);
        var h_fd = filter.GetTransmissionCoefficient(fd / 2);

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

        h_f0_db.AssertLessOrEqualsThan(-Rs);
        h_sl_db.AssertLessOrEqualsThan(-Rp);

        h_pl_db.AssertGreaterOrEqualsThan(-Rs);
        h_c0_db.AssertGreaterOrEqualsThan(-Rp);
        h_ph_db.AssertGreaterOrEqualsThan(-Rs);

        h_sh_db.AssertLessOrEqualsThan(-Rp);
        h_fd_db.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
    public void ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.00000689205724576, +0.00008378176661965, +0.00044120384234372, +0.00123307356096349, +0.00150015430039389,
            -0.00172925558611863, -0.01042871843576654, -0.01781553307525705, -0.00664110958605211, +0.03239945450980376,
            +0.07329424593775179, +0.06015895987979439, -0.03418410317947584, -0.15062005298434142, -0.17296289658174122,
            -0.04432462410435293, +0.14993025382900113, +0.24494314675091355, +0.15750002893034357, -0.02845591991255250,
            -0.15124385623069050, -0.13597042929461414, -0.05084193439842321, -0.00522558571215526, -0.02298973430596381,
            -0.03721003823306662, +0.00706566192563501, +0.08029122954836282, +0.10965905270856308, +0.06984645196828199,
            +0.00750822370827765, -0.02137442334708500, -0.01457963379716382, -0.01221990233503440, -0.03664172372781226,
            -0.06458255123674041, -0.06119874583984340, -0.02539319230568635, +0.01105866175263433, +0.02353771435661804,
            +0.01983529262616530, +0.02248384511405094, +0.03655890545440398, +0.04489236278662814, +0.03269154459100252,
            +0.00687967196457401, -0.01324932462631035, -0.01891056992161094, -0.01859659750849560, -0.02309643407425087,
            -0.03060957592715382, -0.03042136828398954, -0.01794144854210539, -0.00078247890196222, +0.01060663077390213,
            +0.01431820466064217, +0.01636393079078733, +0.02078201057570292, +0.02417417774592334, +0.02084997237448213,
            +0.01069791641060571, -0.00044388844379984, -0.00750085357861182, -0.01090820419717378, -0.01400612559620623,
            -0.01754999037194923, -0.01862081351682704, -0.01471617256857150, -0.00720096643241412, +0.00013320652211490,
            +0.00507824922316659, +0.00844303884552064, +0.01168151220389813, +0.01425566568132750, +0.01420941741997110,
            +0.01077412207082848, +0.00544808428796851, +0.00042538510032214, -0.00341064021649575, -0.00658926881200260,
            -0.00948586329074996, -0.01129283696183008, -0.01085876276739285, -0.00816773063743472, -0.00443734353333156,
            -0.00085243485537876, +0.00228770274777045, +0.00514501475245081, +0.00753411402669440, +0.00878687487962500,
            +0.00834731934573546, +0.00639071248255427, +0.00372761256569559, +0.00105020368895634, -0.00151069458777567,
            -0.00395861646463808, -0.00590495083955218, -0.00678901689076000, -0.00642692800854867, -0.00509018551350081
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

        var filter = new DSP.Filters.ButterworthBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

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

        var filter = new DSP.Filters.ButterworthBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

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
        k_pl_sl_db.AssertGreaterOrEqualsThan(-Rs);         // Коэффициент передачи у верхнего края переходной полосы должен быть гарантировано меньше коэффициента пропускания Rp (-1дБ) и должен приближаться к уровню Rs (-40дБ)

        // Коэффициенты передачи в полосе заграждения
        // должны бытьниже уровня подавления Rs
        k_pl_db.AssertGreaterOrEqualsThan(-Rs);      // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_pl_ph_db.AssertGreaterOrEqualsThan(-Rp);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertGreaterOrEqualsThan(-Rp);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_ph_pl_db.AssertGreaterOrEqualsThan(-Rp);    // Коэффициент передачи у верхнего края полосы подавления
        k_ph_db.AssertGreaterOrEqualsThan(-Rp, 5e-3);      // Коэффициент передачи на верхней границе полосы подавления

        // Коэффициенты передачи в переходной полосе
        // между верхней полосой заграждения и верхней
        // полосой пропускания
        k_ph_sh_db.AssertLessThan(-Rp);
        k_sh_ph_db.AssertLessOrEqualsThan(-Rp);
        k_sh_ph_db.AssertGreaterOrEqualsThan(-Rs);

        // Коэффициенты передачи в верхней полосе пропускания
        // от верхней границы полосы пропускания до частоты,
        // близкой к половине частоты дискретизации
        k_sh_db.AssertLessOrEqualsThan(-Rs, 2.2e-1);
        k_sh_fd05_db.AssertLessOrEqualsThan(-Rs);
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
        h_pl_sl.AssertGreaterOrEqualsThan(-Rs);

        h_pl.AssertGreaterOrEqualsThan(-Rp);
        h_pl_ph.AssertGreaterOrEqualsThan(-Rp);
        h_c0.AssertGreaterOrEqualsThan(-Rp);
        h_ph_pl.AssertGreaterOrEqualsThan(-Rp);
        h_ph.AssertGreaterOrEqualsThan(-Rp, 4e-5);

        h_ph_sh.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertGreaterOrEqualsThan(-Rs);

        h_sh.AssertLessOrEqualsThan(-Rs);
        h_sh_fd05.AssertLessOrEqualsThan(-Rs);
        h_fd05.AssertLessOrEqualsThan(-Rs);
    }
}