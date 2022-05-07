using System;
using System.Diagnostics;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;
using static MathCore.Polynom.Array;
// ReSharper disable InconsistentNaming
// ReSharper disable HeuristicUnreachableCode

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ChebyshevBandPass : ChebyshevFiltersTests
{
    [TestMethod]
    public void TypeI_Even_Creation()
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

        //(Wsl, Wpl, Wph, Wsh).ToDebug();
        //Wsl.ToDebug();
        //Wpl.ToDebug();
        //Wph.ToDebug();
        //Wsh.ToDebug();

        var Wc = Wpl * Wph;
        var dW = Wph - Wpl;

        Wc.AssertEquals(55.472558684693745);

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
        const double F0 = 1 / Consts.pi2;
        var W1 = Abs((Wc - Wp.Pow2()) / (dW * Wp));   // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F1 = W1 / Consts.pi2;

        W1.AssertEquals(1.6258593550796006);
        F1.AssertEquals(0.2587635531331195);

        var eps_p = (1 / (Gp * Gp) - 1).Sqrt();
        //var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.5088471399095873);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.5128464567198);
        kW.AssertEquals(1.6258593550796006);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(6);

        var beta = arcsh(1 / eps_p) / N;
        beta.AssertEquals(0.23799589314392092);

        //var sh = Sinh(beta);
        //var ch = Cosh(beta);

        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -Sinh(beta);
        var r = N % 2;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (im, re) = Complex.Trigonometry.Cos(new(dth * (i - r + 1), -beta));
            (poles[i], poles[i + 1]) = Complex.Conjugate(-re, im);
            //var (cos, sin) = Complex.Exp(dth * (i - r + 1));
            //(poles[i], poles[i + 1]) = Complex.Conjugate(-sh * sin, ch * cos);
        }

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.062181023793011402, 0.993411202482325950),
            /*[ 1]*/ (-0.062181023793011402, -0.993411202482325950),
            /*[ 2]*/ (-0.169881716269156408, 0.727227473025156224),
            /*[ 3]*/ (-0.169881716269156408, -0.727227473025156224),
            /*[ 4]*/ (-0.232062740062167788, 0.266183729457169949),
            /*[ 5]*/ (-0.232062740062167788, -0.266183729457169949)
        );

        var ppf_poles = AnalogBasedFilter.TransformToBandPassW(poles, Wpl, Wph).ToArray();

        //ppf_poles.ToDebugEnum();
        ppf_poles.AssertEquals(
            /*[ 0]*/ (-0.137541851731607462, -4.065157801837609775),
            /*[ 1]*/ (-0.461170340590846606, 13.630252787982325202),
            /*[ 2]*/ (-0.137541851731607462, 4.065157801837609775),
            /*[ 3]*/ (-0.461170340590846606, -13.630252787982325202),
            /*[ 4]*/ (-0.468517038482667514, -4.695467081575139368),
            /*[ 5]*/ (-1.167195090053226414, 11.697602590655399979),
            /*[ 6]*/ (-0.468517038482667514, 4.695467081575139368),
            /*[ 7]*/ (-1.167195090053226414, -11.697602590655399979),
            /*[ 8]*/ (-0.925729759689586484, -6.195367959752425158),
            /*[ 9]*/ (-1.308694561168761705, 8.758327436816884415),
            /*[10]*/ (-0.925729759689586484, 6.195367959752425158),
            /*[11]*/ (-1.308694561168761705, -8.758327436816884415)
        );

        //var z_zeros_enum = DigitalFilter.ToZ(ppf_zeros, dt);
        //if (N.IsOdd())
        //    z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = Enumerable
           .Repeat(Complex.Real, N)
           .Concat(Enumerable.Repeat(-Complex.Real, N))
           .ToArray();
        var z_poles = DigitalFilter.ToZArray(ppf_poles, dt);

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.908563245771859984, -0.385280925843629329),
            /*[ 1]*/ (0.354050639576295845, 0.902003756282745983),
            /*[ 2]*/ (0.908563245771859984, 0.385280925843629329),
            /*[ 3]*/ (0.354050639576295845, -0.902003756282745983),
            /*[ 4]*/ (0.856522835267100269, -0.425885365451733000),
            /*[ 5]*/ (0.447616493057107101, 0.799994631665565725),
            /*[ 6]*/ (0.856522835267100269, 0.425885365451733000),
            /*[ 7]*/ (0.447616493057107101, -0.799994631665565725),
            /*[ 8]*/ (0.757472414797518478, -0.520325379990321801),
            /*[ 9]*/ (0.605874047674048288, 0.660048446020072066),
            /*[10]*/ (0.757472414797518478, 0.520325379990321801),
            /*[11]*/ (0.605874047674048288, -0.660048446020072066)
        );

        var Fp0 = (Fpl * Fph).Sqrt().AssertEquals(1.1853844635393842);
        //var Wp0 = (Wpl * Wph).Sqrt();
        var ffp0 = DigitalFilter.ToDigitalFrequency(Fp0, dt).AssertEquals(1.1347392325852204);
        //var ffp0 = DigitalFilter.ToDigitalFrequency(Wp0 / Consts.pi2, dt).AssertEquals(1.1347392325852204);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);
        z0.AssertEquals(new Complex(0.7564175596225313, 0.6540890424817513));

        var norm_0 = ((z0 - 1) * (z0 + 1)).Pow(N).Abs;
        var norm_p = z_poles.Multiply(z => z0 - z);

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * (norm_p / norm_0).Abs;
        else
            g_norm = (z0 * norm_p / norm_0).Abs;

        g_norm.AssertEquals(0.00018828482383731707);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        filter.B.AssertEquals(Accuracy.Eps(1e-17), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-13), A);
    }

    [TestMethod]
    public void TypeI_Odd_Creation()
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
        const double fph = 12.5 / Consts.pi2; // верхняя частота границы полосы заграждения
        const double fsh = 15 / Consts.pi2; // верхняя частота границы полосы пропускания

        //const double wsl = fsl * Consts.pi2; //  2
        //const double wpl = fpl * Consts.pi2; //  4
        //const double wph = fph * Consts.pi2; // 12
        //const double wsh = fsh * Consts.pi2; // 15

        //const double wc = wpl * wph; // 48
        //const double dw = wph - wpl; // 8

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fsl = DigitalFilter.ToAnalogFrequency(fsl, dt);
        var Fpl = DigitalFilter.ToAnalogFrequency(fpl, dt);
        var Fph = DigitalFilter.ToAnalogFrequency(fph, dt);
        var Fsh = DigitalFilter.ToAnalogFrequency(fsh, dt);

        Fsl.AssertEquals(0.31937518051807723);
        Fpl.AssertEquals(0.64524608331077715);
        Fph.AssertEquals(2.296556302951906);
        Fsh.AssertEquals(2.9653636313402);

        var Wsl = Consts.pi2 * Fsl;
        var Wpl = Consts.pi2 * Fpl;
        var Wph = Consts.pi2 * Fph;
        var Wsh = Consts.pi2 * Fsh;

        //(Wsl, Wpl, Wph, Wsh).ToDebug();
        //Wsl.ToDebug();
        //Wpl.ToDebug();
        //Wph.ToDebug();
        //Wsh.ToDebug();

        var Wc = Wpl * Wph;
        var dW = Wph - Wpl;

        Wc.AssertEquals(58.500854660888386);

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
        const double F0 = 1 / Consts.pi2;
        var W1 = Abs((Wc - Wp.Pow2()) / (dW * Wp));   // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F1 = W1 / Consts.pi2;

        W1.AssertEquals(1.4931453518033904);
        F1.AssertEquals(0.23764146349419665);

        var eps_p = (1 / (Gp * Gp) - 1).Sqrt();
        //var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.5088471399095873);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.5128464567198);
        kW.AssertEquals(1.4931453518033904);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(7);

        var beta = arcsh(1 / eps_p) / N;
        beta.AssertEquals(0.20399647983764652);

        //var sh = Sinh(beta);
        //var ch = Cosh(beta);

        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -Sinh(beta);
        var r = N % 2;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (im, re) = Complex.Trigonometry.Cos(new(dth * (i - r + 1), -beta));
            (poles[i], poles[i + 1]) = Complex.Conjugate(-re, im);
            //var (cos, sin) = Complex.Exp(dth * (i - r + 1));
            //(poles[i], poles[i + 1]) = Complex.Conjugate(-sh * sin, ch * cos);
        }

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/  -0.205414297471464835,
            /*[ 1]*/ (-0.045708981321330548, 0.995283957764552474),
            /*[ 2]*/ (-0.045708981321330548, -0.995283957764552474),
            /*[ 3]*/ (-0.128073719629434535, 0.798155763572583443),
            /*[ 4]*/ (-0.128073719629434535, -0.798155763572583443),
            /*[ 5]*/ (-0.185071887043836425, 0.442943031667010434),
            /*[ 6]*/ (-0.185071887043836425, -0.442943031667010434)
        );

        var ppf_poles = AnalogBasedFilter.TransformToBandPassW(poles, Wpl, Wph).ToArray();

        //ppf_poles.ToDebugEnum();
        ppf_poles.AssertEquals(
            /*[ 0]*/ (-1.065636800483094282, 7.573986590319825574),
            /*[ 1]*/ (-1.065636800483095170, -7.573986590319825574),
            /*[ 2]*/ (-0.104422031157416578, -4.062862706349046782),
            /*[ 3]*/ (-0.369830961046017315, 14.389419575855214362),
            /*[ 4]*/ (-0.104422031157416578, 4.062862706349046782),
            /*[ 5]*/ (-0.369830961046017315, -14.389419575855214362),
            /*[ 6]*/ (-0.347388553297580338, -4.537200988155375647),
            /*[ 7]*/ (-0.981438801875578526, 12.818456622747049778),
            /*[ 8]*/ (-0.347388553297580338, 4.537200988155375647),
            /*[ 9]*/ (-0.981438801875578526, -12.818456622747049778),
            /*[10]*/ (-0.682010245260910919, -5.635385210813174339),
            /*[11]*/ (-1.238200918191908917, 10.231135369124187307),
            /*[12]*/ (-0.682010245260910919, 5.635385210813174339),
            /*[13]*/ (-1.238200918191908917, -10.231135369124187307)
        );

        //var z_zeros_enum = DigitalFilter.ToZ(ppf_zeros, dt);
        //if (N.IsOdd())
        //    z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = Enumerable
           .Repeat(Complex.Real, N)
           .Concat(Enumerable.Repeat(-Complex.Real, N))
           .ToArray();
        var z_poles = DigitalFilter.ToZArray(ppf_poles, dt);

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.681463190332515234, 0.604557069710933370),
            /*[ 1]*/ (0.681463190332515234, -0.604557069710933370),
            /*[ 2]*/ (0.911545352623017502, -0.386300402599526305),
            /*[ 3]*/ (0.309987704704586919, 0.925386310679421231),
            /*[ 4]*/ (0.911545352623017502, 0.386300402599526305),
            /*[ 5]*/ (0.309987704704586919, -0.925386310679421231),
            /*[ 6]*/ (0.872735899271796289, -0.417595562716719493),
            /*[ 7]*/ (0.388272787461148294, 0.848155108648701472),
            /*[ 8]*/ (0.872735899271796289, 0.417595562716719493),
            /*[ 9]*/ (0.388272787461148294, -0.848155108648701472),
            /*[10]*/ (0.800380419932356091, -0.490563396498134574),
            /*[11]*/ (0.528650358164437839, 0.736400828237951277),
            /*[12]*/ (0.800380419932356091, 0.490563396498134574),
            /*[13]*/ (0.528650358164437839, -0.736400828237951277)
        );

        var Fp0 = (Fpl * Fph).Sqrt().AssertEquals(1.2173101328677076);
        //var Wp0 = (Wpl * Wph).Sqrt();
        var ffp0 = DigitalFilter.ToDigitalFrequency(Fp0, dt).AssertEquals(1.1626842508705897);
        //var ffp0 = DigitalFilter.ToDigitalFrequency(Wp0 / Consts.pi2, dt).AssertEquals(1.1626842508705897);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);
        z0.AssertEquals(new Complex(0.7448168130279443, 0.6672689975046767));

        var norm_0 = ((z0 - 1) * (z0 + 1)).Pow(N).Abs;
        var norm_p = z_poles.Multiply(z => z0 - z);

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * (norm_p / norm_0).Abs;
        else
            g_norm = (z0 * norm_p / norm_0).Abs;

        g_norm.AssertEquals(5.782790110607029E-05);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
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

        h_f00.Abs.AssertLessOrEqualsThan(Gs);
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp, 6e-9);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs);

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs);

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);
    }

    [TestMethod]
    public void TypeII_Even_Creation()
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

        //(Wsl, Wpl, Wph, Wsh).ToDebug();
        //Wsl.ToDebug();
        //Wpl.ToDebug();
        //Wph.ToDebug();
        //Wsh.ToDebug();

        var Wc = Wpl * Wph;
        var dW = Wph - Wpl;

        Wc.AssertEquals(55.472558684693745);

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
        const double F0 = 1 / Consts.pi2;
        var W1 = Abs((Wc - Wp.Pow2()) / (dW * Wp));   // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F1 = W1 / Consts.pi2;

        W1.AssertEquals(1.6258593550796006);
        F1.AssertEquals(0.2587635531331195);

        var eps_p = (1 / (Gp * Gp) - 1).Sqrt();
        //var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.5088471399095873);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.5128464567198);
        kW.AssertEquals(1.6258593550796006);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(6);

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(0.8830487276017475);

        //var sh = Sinh(beta);
        //var ch = Cosh(beta);

        var r = N % 2;
        var zeros = new Complex[N - r];
        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -1 / Sinh(beta);
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var th = dth * (i - r + 1);
            var cos = Cos(th);
            (poles[i], poles[i + 1]) = (Complex.i / Complex.Trigonometry.Cos(new(th, beta))).Conjugate();
            (zeros[i], zeros[i + 1]) = Complex.ImValue(1 / cos).Conjugate();
        }

        zeros.Sum(z => z.Re * z.Re).AssertEquals(0);
        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ +1.035276180410083,
            /*[ 1]*/ -1.035276180410083,
            /*[ 2]*/ +1.414213562373095,
            /*[ 3]*/ -1.414213562373095,
            /*[ 4]*/ +3.863703305156270,
            /*[ 5]*/ -3.863703305156270
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.133882765352424271, +0.705787088276199515),
            /*[ 1]*/ (-0.133882765352424271, -0.705787088276199515),
            /*[ 2]*/ (-0.471031461322897138, +0.665351902557851704),
            /*[ 3]*/ (-0.471031461322897138, -0.665351902557851704),
            /*[ 4]*/ (-0.903410457807307554, +0.341931454351975084),
            /*[ 5]*/ (-0.903410457807307554, -0.341931454351975084)
        );

        var ppf_zeros = AnalogBasedFilter.TransformToBandPassW(zeros, Wpl, Wph).ToArray();
        var ppf_poles = AnalogBasedFilter.TransformToBandPassW(poles, Wpl, Wph).ToArray();

        ppf_zeros.ToRe().Average(v => v * v).AssertLessThan(1e-30);

        //ppf_zeros.ToIm().ToDebugEnum();
        ppf_zeros.ToIm().AssertEquals(
            /*[ 0]*/ +13.9458922090638920,
            /*[ 1]*/ -03.9776987985494614,
            /*[ 2]*/ +03.9776987985494614,
            /*[ 3]*/ -13.9458922090638920,
            /*[ 4]*/ +16.8993345941994240,
            /*[ 5]*/ -03.2825291656000655,
            /*[ 6]*/ +03.2825291656000655,
            /*[ 7]*/ -16.8993345941994240,
            /*[ 8]*/ +38.6375215812085300,
            /*[ 9]*/ -01.4357173134954130,
            /*[10]*/ +01.4357173134954130,
            /*[11]*/ -38.6375215812085300
        );
        
        //ppf_poles.ToDebugEnum();
        ppf_poles.AssertEquals(
            /*[ 0]*/ (-0.376334851871586717, -04.767592971414819303),
            /*[ 1]*/ (-0.912760101359974430, +11.563288975735309450),
            /*[ 2]*/ (-0.376334851871586717, +04.767592971414819303),
            /*[ 3]*/ (-0.912760101359974430, -11.563288975735309450),
            /*[ 4]*/ (-1.341047522629964295, -4.6357739010854572830),
            /*[ 5]*/ (-3.194295603920112026, +11.042138286020186655),
            /*[ 6]*/ (-1.341047522629964295, +04.635773901085457283),
            /*[ 7]*/ (-3.194295603920112026, -11.042138286020186655),
            /*[ 8]*/ (-3.224678049868438201, -4.7202458272300038540),
            /*[ 9]*/ (-5.473841575047024399, +08.012544959205451178),
            /*[10]*/ (-3.224678049868438201, +04.720245827230003854),
            /*[11]*/ (-5.473841575047024399, -08.012544959205451178)
        );

        //var z_zeros_enum = DigitalFilter.ToZ(ppf_zeros, dt);
        //if (N.IsOdd())
        //    z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = DigitalFilter.ToZArray(ppf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(ppf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (+0.345695996851573095, +0.938346565913040997),
            /*[ 1]*/ (+0.923899724431281544, -0.382634681119997433),
            /*[ 2]*/ (+0.923899724431281544, +0.382634681119997433),
            /*[ 3]*/ (+0.345695996851573095, -0.938346565913040997),
            /*[ 4]*/ (+0.166882401808782443, +0.985976807012483469),
            /*[ 5]*/ (+0.947538200802179098, -0.319642547262671495),
            /*[ 6]*/ (+0.947538200802179098, +0.319642547262671495),
            /*[ 7]*/ (+0.166882401808782443, -0.985976807012483469),
            /*[ 8]*/ (-0.577358698459896402, +0.816490620468290684),
            /*[ 9]*/ (+0.989746417845975723, -0.142835669078347299),
            /*[10]*/ (+0.989746417845975723, +0.142835669078347299),
            /*[11]*/ (-0.577358698459896402, -0.816490620468290684)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.861171341453314532, -0.435471220443577001),
            /*[ 1]*/ (0.464855504018181320, +0.809962310979558175),
            /*[ 2]*/ (0.861171341453314532, +0.435471220443577001),
            /*[ 3]*/ (0.464855504018181320, -0.809962310979558175),
            /*[ 4]*/ (0.789865665496933267, -0.388800621420344661),
            /*[ 5]*/ (0.405919027746276584, +0.669317688643135011),
            /*[ 6]*/ (0.789865665496933267, +0.388800621420344661),
            /*[ 7]*/ (0.405919027746276584, -0.669317688643135011),
            /*[ 8]*/ (0.653983817737938011, -0.336160105092508654),
            /*[ 9]*/ (0.428872085879188181, +0.449437584642729093),
            /*[10]*/ (0.653983817737938011, +0.336160105092508654),
            /*[11]*/ (0.428872085879188181, -0.449437584642729093)
        );

        var Fp0 = (Fpl * Fph).Sqrt().AssertEquals(1.1853844635393842);
        //var Wp0 = (Wpl * Wph).Sqrt();
        var ffp0 = DigitalFilter.ToDigitalFrequency(Fp0, dt).AssertEquals(1.1347392325852204);
        //var ffp0 = DigitalFilter.ToDigitalFrequency(Wp0 / Consts.pi2, dt).AssertEquals(1.1347392325852204);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);
        z0.AssertEquals(new Complex(0.7564175596225313, 0.6540890424817513));

        var norm_0 = z_zeros.Multiply(z => z0 - z);
        var norm_p = z_poles.Multiply(z => z0 - z);

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * (norm_p / norm_0).Re;
        else
            g_norm = (z0 * norm_p / norm_0).Abs;

        g_norm.AssertEquals(0.014978759442649942);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
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

        h_f00.Abs.AssertLessOrEqualsThan(Gs);
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fpl.Abs.AssertLessOrEqualsThan(Gp);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp, 3.4e-12);
        h_fph.Abs.AssertLessOrEqualsThan(Gp);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs);

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        filter.B.AssertEquals(Accuracy.Eps(1e-15), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-15), A);
    }

    [TestMethod]
    public void TypeII_Odd_Creation()
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
        const double fph = 12.5 / Consts.pi2; // верхняя частота границы полосы заграждения
        const double fsh = 15 / Consts.pi2; // верхняя частота границы полосы пропускания

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fsl = DigitalFilter.ToAnalogFrequency(fsl, dt);
        var Fpl = DigitalFilter.ToAnalogFrequency(fpl, dt);
        var Fph = DigitalFilter.ToAnalogFrequency(fph, dt);
        var Fsh = DigitalFilter.ToAnalogFrequency(fsh, dt);

        Fsl.AssertEquals(0.31937518051807723);
        Fpl.AssertEquals(0.64524608331077715);
        Fph.AssertEquals(2.296556302951906);
        Fsh.AssertEquals(2.9653636313402);

        var Wsl = Consts.pi2 * Fsl;
        var Wpl = Consts.pi2 * Fpl;
        var Wph = Consts.pi2 * Fph;
        var Wsh = Consts.pi2 * Fsh;

        //(Wsl, Wpl, Wph, Wsh).ToDebug();
        //Wpl.ToDebug();
        //Wph.ToDebug();
        //Wsh.ToDebug();

        var Wc = Wpl * Wph;
        var dW = Wph - Wpl;

        Wc.AssertEquals(58.500854660888386);

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
        const double F0 = 1 / Consts.pi2;
        var W1 = Abs((Wc - Wp.Pow2()) / (dW * Wp));   // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F1 = W1 / Consts.pi2;

        W1.AssertEquals(1.4931453518033904);
        F1.AssertEquals(0.23764146349419665);

        var eps_p = (1 / (Gp * Gp) - 1).Sqrt();
        //var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.5088471399095873);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.5128464567198);
        kW.AssertEquals(1.4931453518033904);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(7);

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(0.7568989093729265);

        var r = N % 2;
        var zeros = new Complex[N - r];
        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -1 / Sinh(beta);
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var th = dth * (i - r + 1);
            var cos = Cos(th);
            (poles[i], poles[i + 1]) = (Complex.i / Complex.Trigonometry.Cos(new(th, beta))).Conjugate();
            (zeros[i - r], zeros[i - r + 1]) = Complex.ImValue(1 / cos).Conjugate();
        }

        zeros.Sum(z => z.Re * z.Re).AssertEquals(0);
        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ +1.0257168632725540,
            /*[ 1]*/ -1.0257168632725540,
            /*[ 2]*/ +1.2790480076899327,
            /*[ 3]*/ -1.2790480076899327,
            /*[ 4]*/ +2.3047648709624860,
            /*[ 5]*/ -2.3047648709624860
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/  -1.202981071917573530,
            /*[ 1]*/ (-0.112686910985844033, +0.772336560883440271),
            /*[ 2]*/ (-0.112686910985844033, -0.772336560883440271),
            /*[ 3]*/ (-0.397988375745020884, +0.780702692783778240),
            /*[ 4]*/ (-0.397988375745020884, -0.780702692783778240),
            /*[ 5]*/ (-0.851790250586702546, +0.641693653548654819),
            /*[ 6]*/ (-0.851790250586702546, -0.641693653548654819)
        );

        var ppf_zeros_enum = AnalogBasedFilter.TransformToBandPassW(zeros, Wpl, Wph);
        if (N.IsOdd())
            ppf_zeros_enum = ppf_zeros_enum.AppendLast(0);
        var ppf_zeros = ppf_zeros_enum.ToArray();
        var ppf_poles = AnalogBasedFilter.TransformToBandPassW(poles, Wpl, Wph).ToArray();

        ppf_zeros.ToRe().Average(v => v * v).AssertLessThan(1e-30);

        //ppf_zeros.ToIm().ToDebugEnum();
        ppf_zeros.ToIm().AssertEquals(
            /*[ 0]*/ +14.638643441051567,
            /*[ 1]*/ -03.996330322305191,
            /*[ 2]*/ +03.996330322305191,
            /*[ 3]*/ -14.638643441051567,
            /*[ 4]*/ +16.761036084422237,
            /*[ 5]*/ -03.490288688970682,
            /*[ 6]*/ +03.490288688970682,
            /*[ 7]*/ -16.761036084422237,
            /*[ 8]*/ +26.150172159366733,
            /*[ 9]*/ -02.2371116451688042,
            /*[10]*/ +02.2371116451688042,
            /*[11]*/ -26.150172159366733,
            /*[12]*/ +00
        );

        //ppf_poles.ToDebugEnum();
        ppf_poles.AssertEquals(
            /*[ 0]*/ (-6.240757903904171400, +04.421967372759097792),
            /*[ 1]*/ (-6.240757903904171400, -04.421967372759097792),
            /*[ 2]*/ (-0.312832901614768955, -04.612276119822324105),
            /*[ 3]*/ (-0.856348803431439154, +12.625644923912288320),
            /*[ 4]*/ (-0.312832901614768955, +04.612276119822324105),
            /*[ 5]*/ (-0.856348803431439154, -12.625644923912288320),
            /*[ 6]*/ (-1.076556261828291339, -4.412630915313351387),
            /*[ 7]*/ (-3.052767398490954243, +12.512802421458992796),
            /*[ 8]*/ (-1.076556261828291339, +04.412630915313351387),
            /*[ 9]*/ (-3.052767398490954243, -12.512802421458992796),
            /*[10]*/ (-2.418165073347068983, -04.023548397347719074),
            /*[11]*/ (-6.419574543526488597, +10.681433269776210437),
            /*[12]*/ (-2.418165073347068983, +04.023548397347719074),
            /*[13]*/ (-6.419574543526488597, -10.681433269776210437)
        );

        var z_zeros_enum = DigitalFilter.ToZ(ppf_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendLast(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = DigitalFilter.ToZArray(ppf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (+0.302316746066740527, +0.953207524649075433),
            /*[ 1]*/ (+0.923212583325833669, -0.384289638149196366),
            /*[ 2]*/ (+0.923212583325833669, +0.384289638149196366),
            /*[ 3]*/ (+0.302316746066740527, -0.953207524649075433),
            /*[ 4]*/ (+0.174859768029673179, +0.984593348304064531),
            /*[ 5]*/ (+0.940889647148106745, -0.338713259099066744),
            /*[ 6]*/ (+0.940889647148106745, +0.338713259099066744),
            /*[ 7]*/ (+0.174859768029673179, -0.984593348304064531),
            /*[ 8]*/ (-0.261877886854459507, +0.965101016669572820),
            /*[ 9]*/ (+0.975285872474875637, -0.220946751392548296),
            /*[10]*/ (+0.975285872474875637, +0.220946751392548296),
            /*[11]*/ (-0.261877886854459507, -0.965101016669572820),
            /*[12]*/  +1,
            /*[13]*/  -1
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.482254078640661166, +0.249782388065586641),
            /*[ 1]*/ (0.482254078640661166, -0.249782388065586641),
            /*[ 2]*/ (0.872649850307064989, -0.425207957314191975),
            /*[ 3]*/ (0.403536506475400014, +0.849647929056077444),
            /*[ 4]*/ (0.872649850307064989, +0.425207957314191975),
            /*[ 5]*/ (0.403536506475400014, -0.849647929056077444),
            /*[ 6]*/ (0.818149282152854318, -0.380651451376455618),
            /*[ 7]*/ (0.340276572049154902, +0.727488185963271250),
            /*[ 8]*/ (0.818149282152854318, +0.380651451376455618),
            /*[ 9]*/ (0.340276572049154902, -0.727488185963271250),
            /*[10]*/ (0.728586054103568581, -0.310241700197617043),
            /*[11]*/ (0.301317550140288659, +0.526122650147597293),
            /*[12]*/ (0.728586054103568581, +0.310241700197617043),
            /*[13]*/ (0.301317550140288659, -0.526122650147597293)
        );

        var Fp0 = (Fpl * Fph).Sqrt().AssertEquals(1.2173101328677076);
        //var Wp0 = (Wpl * Wph).Sqrt();
        var ffp0 = DigitalFilter.ToDigitalFrequency(Fp0, dt).AssertEquals(1.1626842508705897);
        //var ffp0 = DigitalFilter.ToDigitalFrequency(Wp0 / Consts.pi2, dt).AssertEquals(1.1347392325852204);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);
        z0.AssertEquals(new Complex(0.7448168130279443, 0.6672689975046767));

        var norm_0 = z_zeros.Multiply(z => z0 - z);
        var norm_p = z_poles.Multiply(z => z0 - z);

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * (norm_p / norm_0).Re;
        else
            g_norm = (norm_p / norm_0).Re;

        g_norm.AssertEquals(0.018631883267691787);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
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

        h_f00.Abs.AssertLessOrEqualsThan(Gs);
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fpl.Abs.AssertLessOrEqualsThan(Gp);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fph.Abs.AssertLessOrEqualsThan(Gp);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs);

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(Accuracy.Eps(1e-13), A);
    }

    [TestMethod, Ignore]
    public void TypeIICorrected_Even_Creation()
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
        //const double f1 = 1 / Consts.pi2;   // 0.159
        //var w1 = f1 * Consts.pi2;

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

        var Kw = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(kW);
        var Tw = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(kW);
        var K_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(kEps);
        var T_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(kEps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N);
        N.AssertEquals(4);

        //var Fp = DigitalFilter.ToAnalogFrequency(f0, dt);
        //var Fs = DigitalFilter.ToAnalogFrequency(f1, dt);

        var (L, r) = N.GetDivMod(2);
        (L, r).AssertEquals((2, 0));

        var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);
        var m = (1 - kEps * kEps).Sqrt();

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * SpecialFunctions.EllipticJacobi.sn_uk(ui, m).Power(4));

        var k_W = Sqrt(1 - kp * kp);

        //var im_pz = Range(0, L).ToArray(i => 1 / (k_W * cd_uk(u[i], k_W)));

        var v0_complex = SpecialFunctions.EllipticJacobi.sn_inverse((0, 1 / eps_p), kEps) / N;

        var zeros = new Complex[N - r];
        var poles = new Complex[N];

        if (r != 0) poles[0] = Complex.i * SpecialFunctions.EllipticJacobi.sn_uk(v0_complex, k_W);
        for (var i = 0; i < L; i++)
        {
            var (p_im, p_re) = SpecialFunctions.EllipticJacobi.cd_uk(u[i] - v0_complex, k_W);
            (poles[r + 2 * i], poles[r + 2 * i + 1]) = Complex.Conjugate(-p_re, p_im);

            var p0_im = 1 / (k_W * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k_W));
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

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * (norm_p / norm_0).Abs;
        else
            g_norm = (z0 * norm_p / norm_0).Abs;

        g_norm.AssertEquals(0.027089200894329788);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
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

    [TestMethod, Ignore]
    public void TypeIICorrected_Odd_Creation()
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
        //const double f1 = 1 / Consts.pi2;   // 0.159
        //var w1 = f1 * Consts.pi2;

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

        var Kw = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(kW);
        var Tw = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(kW);
        var K_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(kEps);
        var T_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(kEps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N);
        N.AssertEquals(4);

        //var Fp = DigitalFilter.ToAnalogFrequency(f0, dt);
        //var Fs = DigitalFilter.ToAnalogFrequency(f1, dt);

        var (L, r) = N.GetDivMod(2);
        (L, r).AssertEquals((2, 0));

        var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);
        var m = (1 - kEps * kEps).Sqrt();

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * SpecialFunctions.EllipticJacobi.sn_uk(ui, m).Power(4));

        var k_W = Sqrt(1 - kp * kp);

        //var im_pz = Range(0, L).ToArray(i => 1 / (k_W * cd_uk(u[i], k_W)));

        var v0_complex = SpecialFunctions.EllipticJacobi.sn_inverse((0, 1 / eps_p), kEps) / N;

        var zeros = new Complex[N - r];
        var poles = new Complex[N];

        if (r != 0) poles[0] = Complex.i * SpecialFunctions.EllipticJacobi.sn_uk(v0_complex, k_W);
        for (var i = 0; i < L; i++)
        {
            var (p_im, p_re) = SpecialFunctions.EllipticJacobi.cd_uk(u[i] - v0_complex, k_W);
            (poles[r + 2 * i], poles[r + 2 * i + 1]) = Complex.Conjugate(-p_re, p_im);

            var p0_im = 1 / (k_W * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k_W));
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

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * (norm_p / norm_0).Abs;
        else
            g_norm = (z0 * norm_p / norm_0).Abs;

        g_norm.AssertEquals(0.027089200894329788);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
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

    [TestMethod, Ignore]
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

        h_f0_db.AssertLessOrEqualsThan(-Rs, 2e-12);
        h_sl_db.AssertLessOrEqualsThan(-Rp);

        h_pl_db.AssertGreaterOrEqualsThan(-Rp);
        h_c0_db.AssertGreaterOrEqualsThan(-Rp);
        h_ph_db.AssertGreaterOrEqualsThan(-Rp);

        h_sh_db.AssertLessOrEqualsThan(-Rp);
        h_fd_db.AssertLessOrEqualsThan(-Rs, 4.1e-12);
    }

    [TestMethod, Ignore]
    public void TypeI_Even_ImpulseResponse()
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

    [TestMethod, Ignore]
    public void TypeI_Odd_ImpulseResponse()
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

    [TestMethod, Ignore]
    public void TypeII_Even_ImpulseResponse()
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

    [TestMethod, Ignore]
    public void TypeII_Odd_ImpulseResponse()
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

    [TestMethod, Ignore]
    public void TypeIICorrected_Even_ImpulseResponse()
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

    [TestMethod, Ignore]
    public void TypeIICorrected_Odd_ImpulseResponse()
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

    [TestMethod, Ignore]
    public void TypeI_Even_SignalProcessing()
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

    [TestMethod, Ignore]
    public void TypeI_Odd_SignalProcessing()
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

    [TestMethod, Ignore]
    public void TypeII_Even_SignalProcessing()
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

    [TestMethod, Ignore]
    public void TypeII_Odd_SignalProcessing()
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

    [TestMethod, Ignore]
    public void TypeIICorrected_Even_SignalProcessing()
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

    [TestMethod, Ignore]
    public void TypeIICorrected_Odd_SignalProcessing()
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