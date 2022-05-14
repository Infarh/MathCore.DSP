using System;
using System.Diagnostics;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Infrastructure;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;
using static MathCore.Polynom.Array;
// ReSharper disable InconsistentNaming
// ReSharper disable HeuristicUnreachableCode
// ReSharper disable RedundantArgumentDefaultValue

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
        poles.AssertEquals(AccuracyComplex.Eps(1e-15),
            /*[ 0]*/ (-0.062181023793011402, 0.993411202482325950),
            /*[ 1]*/ (-0.062181023793011402, -0.993411202482325950),
            /*[ 2]*/ (-0.169881716269156408, 0.727227473025156224),
            /*[ 3]*/ (-0.169881716269156408, -0.727227473025156224),
            /*[ 4]*/ (-0.232062740062167788, 0.266183729457169949),
            /*[ 5]*/ (-0.232062740062167788, -0.266183729457169949)
        );

        var ppf_poles = AnalogBasedFilter.TransformToBandPassW(poles, Wpl, Wph).ToArray();

        //ppf_poles.ToDebugEnum();
        ppf_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
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
        Fpl.AssertEquals(0.64524608331077715, 1e-15);
        Fph.AssertEquals(2.296556302951906, 1e-14);
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

        Wc.AssertEquals(58.500854660888386, 1e-13);

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

        W1.AssertEquals(1.4931453518033904, 1e-15);
        F1.AssertEquals(0.23764146349419665, 1e-15);

        var eps_p = (1 / (Gp * Gp) - 1).Sqrt();
        //var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.5088471399095873);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.5128464567198);
        kW.AssertEquals(1.4931453518033904, 1e-14);

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
        ppf_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
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
        z_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
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
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp, 1.57e-11);
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
        poles.AssertEquals(AccuracyComplex.Eps(1e-15),
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

        h_fcc.Abs.ToDebug();
        h_fcc.Abs.In_dB().ToDebug();

        h_f00.Abs.AssertLessOrEqualsThan(Gs);
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fpl.Abs.AssertLessOrEqualsThan(Gp);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp, 3.4e-12);
        h_fph.Abs.AssertLessOrEqualsThan(Gp);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs);

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        filter.B.AssertEquals(Accuracy.Eps(1e-14), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-10), A);
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
        Fph.AssertEquals(2.296556302951906, 1e-15);
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

        Wc.AssertEquals(58.500854660888386, 1e-12);

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

        W1.AssertEquals(1.4931453518033904, 9e-16);
        F1.AssertEquals(0.23764146349419665, 1e-14);

        var eps_p = (1 / (Gp * Gp) - 1).Sqrt();
        //var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.5088471399095873);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.5128464567198);
        kW.AssertEquals(1.4931453518033904, 9e-16);

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
        ppf_zeros.ToIm().AssertEquals(Accuracy.Eps(1e-14),
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
        ppf_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
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
        z_zeros.AssertEquals(AccuracyComplex.Eps(1e-14),
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
        z_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
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

        g_norm.AssertEquals(0.018631883267691787, 1e-17);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        B.AssertEquals(Accuracy.Eps(1e-14),
            /*[ 0]*/ 0.018631883267691787,
            /*[ 1]*/ -0.11382913315256311,
            /*[ 2]*/ 0.3335771258287947,
            /*[ 3]*/ -0.6352334627094309,
            /*[ 4]*/ 0.8738501682103076,
            /*[ 5]*/ -0.8703867503348532,
            /*[ 6]*/ 0.5449200052977058,
            /*[ 7]*/ -5.295499723755423E-16,
            /*[ 8]*/ -0.5449200052977058,
            /*[ 9]*/ 0.8703867503348532,
            /*[10]*/ -0.8738501682103078,
            /*[11]*/ 0.635233462709431,
            /*[12]*/ -0.3335771258287948,
            /*[13]*/ 0.11382913315256316,
            /*[14]*/ -0.018631883267691784
        );

        //A.ToDebugEnum();
        A.AssertEquals(Accuracy.Eps(1e-10),
            /*[ 0]*/ 1,
            /*[ 1]*/ -7.893539787737986,
            /*[ 2]*/ 30.602554618398255,
            /*[ 3]*/ -76.98535753402601,
            /*[ 4]*/ 140.06725917701425,
            /*[ 5]*/ -194.55773046868126,
            /*[ 6]*/ 212.40595930060715,
            /*[ 7]*/ -184.93636875723766,
            /*[ 8]*/ 128.96380351471478,
            /*[ 9]*/ -71.67444915764588,
            /*[10]*/ 31.27025009748173,
            /*[11]*/ -10.399248150520812,
            /*[12]*/ 2.497555201285449,
            /*[13]*/ -0.3888957362494067,
            /*[14]*/ 0.02977300184748293
        );

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

        h_fcc.Abs.ToDebug();
        h_fcc.Abs.In_dB().ToDebug();

        h_f00.Abs.AssertLessOrEqualsThan(Gs);
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fpl.Abs.AssertLessOrEqualsThan(Gp);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fph.Abs.AssertLessOrEqualsThan(Gp);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs);

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        filter.B.AssertEquals(Accuracy.Eps(1e-10), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-10), A);
    }

    [TestMethod]
    public void TypeII_Corrected_Even_Creation()
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
        var Ws = Wc / Wsh > Wsl
            ? Wsh
            : Wsl;
        const double W0 = 1;
        const double F0 = W0 / Consts.pi2;
        var W1 = Abs((Wc - Ws.Pow2()) / (dW * Ws));   // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
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

        var r = N % 2;
        var sh = Sinh(beta);
        var ch = Cosh(beta);

        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -kW / sh;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (sin, cos) = Complex.SinCos(dth * (i - r + 1));
            var z = new Complex(-sin * sh, cos * ch);
            var norm = kW / z.Power;
            (poles[i], poles[i + 1]) = Complex.Conjugate(z.Re * norm, z.Im * norm);
        }
        //poles.ToDebugEnum();

        var zeros = new Complex[N - r];
        for (var (n, dth, L) = (1, PI / N, N / 2); n <= L; n++)
        {
            var th = dth * (n - 0.5);
            (zeros[2 * n - 2], zeros[2 * n - 1]) = Complex.Conjugate(0, kW / Cos(th));
        }

        zeros.ToRe().Sum(v => v * v).AssertEquals(0);
        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ +1.6832134630108098,
            /*[ 1]*/ -1.6832134630108098,
            /*[ 2]*/ +2.2993123504647450,
            /*[ 3]*/ -2.2993123504647450,
            /*[ 4]*/ +6.2818381639402950,
            /*[ 5]*/ -6.2818381639402950
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(AccuracyComplex.Eps(1e-15),
            /*[ 0]*/ (-0.217674546532166030, +1.147510540168250825),
            /*[ 1]*/ (-0.217674546532166030, -1.147510540168250825),
            /*[ 2]*/ (-0.765830907928647497, +1.081768615193693828),
            /*[ 3]*/ (-0.765830907928647497, -1.081768615193693828),
            /*[ 4]*/ (-1.468818344302755596, +0.555932453854131925),
            /*[ 5]*/ (-1.468818344302755596, -0.555932453854131925)
        );

        var ppf_zeros = AnalogBasedFilter.TransformToBandPassW(zeros, Wpl, Wph).ToArray();
        var ppf_poles = AnalogBasedFilter.TransformToBandPassW(poles, Wpl, Wph).ToArray();

        ppf_zeros.ToRe().Sum(v => v * v).AssertEquals(0, 1e-25);
        //ppf_zeros.ToIm().ToDebugEnum();
        ppf_zeros.ToIm().AssertEquals(
            /*[ 0]*/ 19.10972494356975200,
            /*[ 1]*/ -02.9028444338420343,
            /*[ 2]*/ +02.9028444338420343,
            /*[ 3]*/ -19.1097249435697520,
            /*[ 4]*/ +24.4114130962799360,
            /*[ 5]*/ -02.2724026038929814,
            /*[ 6]*/ +02.2724026038929814,
            /*[ 7]*/ -24.4114130962799360,
            /*[ 8]*/ +61.3885321361604350,
            /*[ 9]*/ -00.9036306416588538,
            /*[10]*/ +00.9036306416588538,
            /*[11]*/ -61.3885321361604350
        );

        //ppf_poles.ToDebugEnum();
        ppf_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/ (-00.421063556293708952, -03.710651515169726800),
            /*[ 1]*/ (-01.674823533003725107, +14.759497438071257136),
            /*[ 2]*/ (-00.421063556293708952, +03.710651515169726800),
            /*[ 3]*/ (-01.674823533003725107, -14.759497438071257136),
            /*[ 4]*/ (-01.454596782213441575, -3.3935255997224302860),
            /*[ 5]*/ (-05.919233268583965923, +13.809373067017329362),
            /*[ 6]*/ (-01.454596782213441575, +03.393525599722430286),
            /*[ 7]*/ (-05.919233268583965923, -13.809373067017329362),
            /*[ 8]*/ (-03.378655139160580134, -02.448839845406638460),
            /*[ 9]*/ (-10.763914368351722217, +07.801655188849368372),
            /*[10]*/ (-03.378655139160580134, +02.448839845406638460),
            /*[11]*/ (-10.763914368351722217, -07.801655188849368372)
        );

        var z_zeros_enum = DigitalFilter.ToZ(ppf_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = DigitalFilter.ToZArray(ppf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (+0.045503463692268660, +0.998964180935435930),
            /*[ 1]*/ (+0.958736733844263278, -0.284295401260087288),
            /*[ 2]*/ (+0.958736733844263278, +0.284295401260087288),
            /*[ 3]*/ (+0.045503463692268660, -0.998964180935435930),
            /*[ 4]*/ (-0.196720280684885118, +0.980459653003253440),
            /*[ 5]*/ (+0.974509996100728748, -0.224344082827600821),
            /*[ 6]*/ (+0.974509996100728748, +0.224344082827600821),
            /*[ 7]*/ (-0.196720280684885118, -0.980459653003253440),
            /*[ 8]*/ (-0.808086831243661985, +0.589063386377542186),
            /*[ 9]*/ (+0.995925575728179302, -0.090178975434928574),
            /*[10]*/ (+0.995925575728180079, +0.090178975434928643),
            /*[11]*/ (-0.808086831243661985, -0.589063386377542186)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.896155602028712472, -0.344544868501569401),
            /*[ 1]*/ (0.260822535846387316, +0.858558638752081493),
            /*[ 2]*/ (0.896155602028712472, +0.344544868501569401),
            /*[ 3]*/ (0.260822535846387316, -0.858558638752081493),
            /*[ 4]*/ (0.818896234695045133, -0.287699228204258417),
            /*[ 5]*/ (0.202043819532995661, +0.640430655290771256),
            /*[ 6]*/ (0.818896234695045133, +0.287699228204258417),
            /*[ 7]*/ (0.202043819532995661, -0.640430655290771256),
            /*[ 8]*/ (0.692393635600003887, -0.177272856128827733),
            /*[ 9]*/ (0.221657713376083088, +0.309809477572961622),
            /*[10]*/ (0.692393635600003887, +0.177272856128827733),
            /*[11]*/ (0.221657713376083088, -0.309809477572961622)
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

        g_norm.AssertEquals(0.027988054007505623);

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
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp, 4.6e-2);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp, 4.6e-2);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs);

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_Creation()
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
        Fph.AssertEquals(2.296556302951906, 1e-15);
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

        Wc.AssertEquals(58.500854660888386, 1e-13);

        // Выбор опорной частоты
        // Если   Wc / Wsh > Wsl
        // то есть      Wc > Wsl*Wsh
        // то есть Wpl*Wph > Wsl*Wsh
        // то есть центральная частота по границам подавления > центральной частоты по границам пропускания
        // то выбираем в качестве опорной частоты выбираем верхнюю границу пропускания
        // иначе, выбираем нижнюю границу пропускания
        var Ws = Wc / Wsh > Wsl
            ? Wsh
            : Wsl;
        const double W0 = 1;
        const double F0 = W0 / Consts.pi2;
        var W1 = Abs((Wc - Ws.Pow2()) / (dW * Ws));   // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F1 = W1 / Consts.pi2;

        W1.AssertEquals(1.4931453518033904, 1e-15);
        F1.AssertEquals(0.23764146349419665, 1e-15);

        var eps_p = (1 / (Gp * Gp) - 1).Sqrt();
        //var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.5088471399095873);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.5128464567198);
        kW.AssertEquals(1.4931453518033904, 1e-13);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(7);

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(0.7568989093729265);

        var r = N % 2;
        var sh = Sinh(beta);
        var ch = Cosh(beta);

        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -kW / sh;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (sin, cos) = Complex.SinCos(dth * (i - r + 1));
            var z = new Complex(-sin * sh, cos * ch);
            var norm = kW / z.Power;
            (poles[i], poles[i + 1]) = Complex.Conjugate(z.Re * norm, z.Im * norm);
        }
        //poles.ToDebugEnum();

        var zeros = new Complex[N - r];
        for (var (n, dth, L) = (1, PI / N, N / 2); n <= L; n++)
        {
            var th = dth * (n - 0.5);
            (zeros[2 * n - 2], zeros[2 * n - 1]) = Complex.Conjugate(0, kW / Cos(th));
        }

        zeros.ToRe().Sum(v => v * v).AssertEquals(0);
        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(Accuracy.Eps(1e-14),
            /*[ 0]*/ +1.5315443666617676,
            /*[ 1]*/ -1.5315443666617676,
            /*[ 2]*/ +1.9098045874156100,
            /*[ 3]*/ -1.9098045874156100,
            /*[ 4]*/ +3.4413489540773770,
            /*[ 5]*/ -3.4413489540773770
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/  -1.796225595841184797,
            /*[ 1]*/ (-0.168257937347595449, +1.153210745910925095),
            /*[ 2]*/ (-0.168257937347595449, -1.153210745910925095),
            /*[ 3]*/ (-0.594254493315459098, +1.165702596870488827),
            /*[ 4]*/ (-0.594254493315459098, -1.165702596870488827),
            /*[ 5]*/ (-1.271846653374979930, +0.958141896077909005),
            /*[ 6]*/ (-1.271846653374979930, -0.958141896077909005)
        );

        var ppf_zeros_enum = AnalogBasedFilter.TransformToBandPassW(zeros, Wpl, Wph);
        if (N.IsOdd())
            ppf_zeros_enum = ppf_zeros_enum.AppendLast(0);
        var ppf_zeros = ppf_zeros_enum.ToArray();
        var ppf_poles = AnalogBasedFilter.TransformToBandPassW(poles, Wpl, Wph).ToArray();

        ppf_zeros.ToRe().Sum(v => v * v).AssertEquals(0, 1e-25);
        //ppf_zeros.ToIm().ToDebugEnum();
        ppf_zeros.ToIm().AssertEquals(Accuracy.Eps(1e-14),
            /*[ 0]*/ +18.9737693311404280,
            /*[ 1]*/ -03.0832489654480346,
            /*[ 2]*/ +03.0832489654480346,
            /*[ 3]*/ -18.9737693311404280,
            /*[ 4]*/ +22.4240043797515440,
            /*[ 5]*/ -02.6088495912761047,
            /*[ 6]*/ +02.6088495912761047,
            /*[ 7]*/ -22.4240043797515440,
            /*[ 8]*/ +37.2751099522368200,
            /*[ 9]*/ -01.5694347980689933,
            /*[10]*/ +01.5694347980689933,
            /*[11]*/ -37.2751099522368200,
            /*[12]*/ +00
        );

        //ppf_poles.ToDebugEnum();
        ppf_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/  -03.995675744110501526,
            /*[ 1]*/  -14.641041567779060628,
            /*[ 2]*/ (-0.3337508338309236300, -03.703543651100382306),
            /*[ 3]*/ (-1.4120073944723849900, +15.668668033213606350),
            /*[ 4]*/ (-0.3337508338309236300, +03.703543651100382306),
            /*[ 5]*/ (-1.4120073944723849900, -15.668668033213606350),
            /*[ 6]*/ (-1.1116896710854820000, -03.410594455978483097),
            /*[ 7]*/ (-5.0539907584119614300, +15.505327889190112955),
            /*[ 8]*/ (-1.1116896710854820000, +03.410594455978483097),
            /*[ 9]*/ (-5.0539907584119614300, -15.505327889190112955),
            /*[10]*/ (-2.2978152938911033900, -02.656041590088078053),
            /*[11]*/ (-10.898214535492323662, +12.597231440196788199),
            /*[12]*/ (-2.2978152938911033940, +02.656041590088078053),
            /*[13]*/ (-10.898214535492323662, -12.597231440196788199)
        );

        var z_zeros_enum = DigitalFilter.ToZ(ppf_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendLast(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = DigitalFilter.ToZArray(ppf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/ (+0.052626145968034464, +0.998614284276241992),
            /*[ 1]*/ (+0.953571306228448101, -0.301167335442891426),
            /*[ 2]*/ (+0.953571306228448101, +0.301167335442891426),
            /*[ 3]*/ (+0.052626145968034464, -0.998614284276241992),
            /*[ 4]*/ (-0.113903273201571217, +0.993491844130573787),
            /*[ 5]*/ (+0.966538868253713157, -0.256520206133613649),
            /*[ 6]*/ (+0.966538868253713157, +0.256520206133613649),
            /*[ 7]*/ (-0.113903273201571217, -0.993491844130573787),
            /*[ 8]*/ (-0.552931217580533541, +0.833226901044912704),
            /*[ 9]*/ (+0.987759745284454427, -0.155982965722509054),
            /*[10]*/ (+0.987759745284454427, +0.155982965722509054),
            /*[11]*/ (-0.552931217580533541, -0.833226901044912704),
            /*[12]*/  +1,
            /*[13]*/  -1
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/  0.666967016330748841,
            /*[ 1]*/  0.154699691166489295,
            /*[ 2]*/ (0.904008828827797872, -0.346791886419301199),
            /*[ 3]*/ (0.216623786479769070, +0.890288980405640307),
            /*[ 4]*/ (0.904008828827797872, +0.346791886419301199),
            /*[ 5]*/ (0.216623786479769070, -0.890288980405640307),
            /*[ 6]*/ (0.846494336983498741, -0.298301246694471001),
            /*[ 7]*/ (0.154405321789832628, +0.714434407036235508),
            /*[ 8]*/ (0.846494336983498741, +0.298301246694471001),
            /*[ 9]*/ (0.154405321789832628, -0.714434407036235508),
            /*[10]*/ (0.768800593285895784, -0.210693643230012412),
            /*[11]*/ (0.110059090550493829, +0.452572146519861163),
            /*[12]*/ (0.768800593285895784, +0.210693643230012412),
            /*[13]*/ (0.110059090550493829, -0.452572146519861163)
        );

        var Fp0 = (Fpl * Fph).Sqrt().AssertEquals(1.2173101328677076);
        var ffp0 = DigitalFilter.ToDigitalFrequency(Fp0, dt).AssertEquals(1.1626842508705897);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);
        z0.AssertEquals(new Complex(0.7448168130279443, 0.6672689975046767));

        var norm_0 = z_zeros.Multiply(z => z0 - z);
        var norm_p = z_poles.Multiply(z => z0 - z);

        //(norm_0 / norm_p).Re.ToDebug();

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * (norm_p / norm_0).Abs;
        else
            g_norm = (z0 * norm_p / norm_0).Abs;

        g_norm.AssertEquals(0.03204624482119541, 1e-16);

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
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp, 4.6e-2);
        h_fcc.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp, 4.6e-2);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs, 5e-14);
        h_fd5.Abs.AssertLessOrEqualsThan(Gs);

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

        filter.B.AssertEquals(Accuracy.Eps(1e-10), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-10), A);
    }

    [TestMethod]
    public void TypeI_Even_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

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

        h_pl_db.AssertGreaterOrEqualsThan(-Rp);
        h_c0_db.AssertGreaterOrEqualsThan(-Rp);
        h_ph_db.AssertGreaterOrEqualsThan(-Rp, 6e-11);

        h_sh_db.AssertLessOrEqualsThan(-Rp);
        h_fd_db.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
    public void TypeI_Odd_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

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

        h_pl_db.AssertGreaterOrEqualsThan(-Rp, 6e-8);
        h_c0_db.AssertGreaterOrEqualsThan(-Rp);
        h_ph_db.AssertGreaterOrEqualsThan(-Rp, 1.6e-10);

        h_sh_db.AssertLessOrEqualsThan(-Rp);
        h_fd_db.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
    public void TypeII_Even_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        var h_f0 = filter.GetTransmissionCoefficient(0);
        var h_sl = filter.GetTransmissionCoefficient(fsl);
        var h_pl = filter.GetTransmissionCoefficient(fpl);
        var h_c0 = filter.GetTransmissionCoefficient((fpl * fph).Sqrt());
        var h_ph = filter.GetTransmissionCoefficient(fph);
        var h_sh = filter.GetTransmissionCoefficient(fsh);
        var h_fd = filter.GetTransmissionCoefficient(fd / 2);

        //fsl.ToDebug();
        //fpl.ToDebug();
        //fph.ToDebug();
        //fsh.ToDebug();

        //h_f0.Abs.In_dB().ToDebug();
        //h_sl.Abs.In_dB().ToDebug();
        //h_pl.Abs.In_dB().ToDebug();
        h_c0.Abs.In_dB().ToDebug();
        //h_ph.Abs.In_dB().ToDebug();
        //h_sh.Abs.In_dB().ToDebug();
        //h_fd.Abs.In_dB().ToDebug();

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
        h_sl_db.AssertLessOrEqualsThan(-Rs);

        h_pl_db.AssertLessOrEqualsThan(-Rs);

        h_c0_db.AssertEquals(-Rp, 4.9e-11);

        h_ph_db.AssertLessOrEqualsThan(-Rs);

        h_sh_db.AssertLessOrEqualsThan(-Rs);
        h_fd_db.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
    public void TypeII_Odd_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        var h_f0 = filter.GetTransmissionCoefficient(0);
        var h_sl = filter.GetTransmissionCoefficient(fsl);
        var h_pl = filter.GetTransmissionCoefficient(fpl);
        var h_c0 = filter.GetTransmissionCoefficient((fpl * fph).Sqrt());
        var h_ph = filter.GetTransmissionCoefficient(fph);
        var h_sh = filter.GetTransmissionCoefficient(fsh);
        var h_fd = filter.GetTransmissionCoefficient(fd / 2);

        //h_f0.Power.In_dB_byPower().ToDebug();
        //h_sl.Power.In_dB_byPower().ToDebug();
        //h_pl.Power.In_dB_byPower().ToDebug();
        //h_c0.Abs.ToDebug();
        //h_ph.Power.In_dB_byPower().ToDebug();
        //h_sh.Power.In_dB_byPower().ToDebug();
        //h_fd.Power.In_dB_byPower().ToDebug();

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
        h_sl_db.AssertLessOrEqualsThan(-Rs);

        h_pl_db.AssertLessOrEqualsThan(-Rs, 6e-9);
        h_c0_db.AssertGreaterOrEqualsThan(-Rp);
        h_c0_db.AssertEquals(0, 1.3e-10);
        h_ph_db.AssertLessOrEqualsThan(-Rs, 2.6e-11);

        h_sh_db.AssertLessOrEqualsThan(-Rs);
        h_fd_db.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
    public void TypeII_Corrected_Even_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

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

        h_pl_db.AssertGreaterOrEqualsThan(-Rp, 0.46);
        h_c0_db.AssertGreaterOrEqualsThan(-Rp, 1e-12);
        h_ph_db.AssertGreaterOrEqualsThan(-Rp, 0.46);

        h_sh_db.AssertLessOrEqualsThan(-Rp);
        h_fd_db.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

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

        h_pl_db.AssertGreaterOrEqualsThan(-Rp);
        h_c0_db.AssertGreaterOrEqualsThan(-Rp);
        h_ph_db.AssertGreaterOrEqualsThan(-Rp);

        h_sh_db.AssertLessOrEqualsThan(-Rp);
        h_fd_db.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
    public void TypeI_Even_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.000188284860342486, +0.001479956305863380, +0.004749902982439790, +0.006626656748170020, -0.002508700952507990,
            -0.027036972447236900, -0.047658930380903000, -0.026875522222909500, +0.048156906403668100, +0.127518952640903000,
            +0.125471055413688000, +0.007199783352445440, -0.153274544241392000, -0.222985771200800000, -0.133943816158011000,
            +0.045014843347039800, +0.170318274178307000, +0.159405560148488000, +0.058976711036568100, -0.020495460062454300,
            -0.024123613028188900, +0.005381232055974190, -0.002406472562804330, -0.056100501168090500, -0.097131880408097200,
            -0.075123193121255600, -0.007892828803182270, +0.042639905274967700, +0.042734268298674500, +0.018681384049698800,
            +0.015751425787796900, +0.040550010009382200, +0.056801780536469400, +0.034122538892406300, -0.013715963681112100,
            -0.045948160711845900, -0.042076589654004700, -0.020758959260145500, -0.011673209094949300, -0.019473763936439400,
            -0.023539106620881500, -0.008750770874050140, +0.014624168091743000, +0.025583569584433500, +0.019550397773783100,
            +0.012234365326072100, +0.016944025107528000, +0.026391792406308700, +0.022572682676533000, +0.001095596925168210,
            -0.021775997147370400, -0.028115701574587000, -0.018918395689351700, -0.010728988422396600, -0.013716200461413900,
            -0.019425093797297700, -0.013011846718691000, +0.006514122496273820, +0.023827839976293800, +0.025283525949237100,
            +0.014432267686419000, +0.006218445501650210, +0.008174708981836740, +0.012176831137656800, +0.006340224134454350,
            -0.008914763914320470, -0.020771674472653600, -0.019371276856041700, -0.009114767408211420, -0.002390764047915550,
            -0.004377544528578730, -0.007877341326967550, -0.003687979784521760, +0.007285462503163000, +0.015172625325745700,
            +0.013317421814009700, +0.005905700026842920, +0.002033447916076970, +0.004384102869773070, +0.006718852935048130,
            +0.002611809965597700, -0.006110425846211280, -0.011644978161238900, -0.009770327310897160, -0.004510704406457110,
            -0.002520749128567200, -0.004776255912809090, -0.005961753341349320, -0.001782503217825400, +0.005365472713167630,
            +0.009318096150187050, +0.007583322840720870, +0.003897148646025110, +0.003080374478408660, +0.005028654108295740,
            +0.005299092144551790, +0.001146184032625430, -0.004811040624158530, -0.007729208292000710, -0.006268888476856120
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod]
    public void TypeI_Odd_ImpulseResponse()
    {
        double[] expected_h =
        {
            +5.7827877333733e-005, +0.000519645569176767, +0.001940177682683280, +0.003362853166507800, +0.000111463882881050,
            -0.012515732982596400, -0.027381321814588900, -0.020802327978941400, +0.025126310539700000, +0.086255814532957000,
            +0.096299248323145800, +0.009212156002258640, -0.130401330094283000, -0.201284751052029000, -0.117511573544446000,
            +0.070891135096577100, +0.209813497749611000, +0.188041104753759000, +0.045103812949553300, -0.079051174566505000,
            -0.094198557190624300, -0.039874456385998300, -0.013837600936063800, -0.047237945512662400, -0.077791963692760100,
            -0.041494082742196500, +0.041533087480998200, +0.092985182877335500, +0.072955759062066600, +0.021921026797681900,
            +0.001494253099410410, +0.017453921488334300, +0.022740623297132400, -0.012025107661415900, -0.057672325144649000,
            -0.065947200445340500, -0.032261968070538700, +0.001244744617363310, +0.002207099318371930, -0.012912636873671800,
            -0.006105465945587130, +0.028808419071681800, +0.057615373634768100, +0.049943767338122200, +0.016057287710447100,
            -0.009450032236347850, -0.011135452252837700, -0.005909641865932880, -0.013696172496154800, -0.029273929112872200,
            -0.032111887189457200, -0.016603545873566100, +0.000416640158841018, +0.003413648982314690, -0.001569922295619740,
            +0.002300443889987560, +0.017514208585507400, +0.028238199404819700, +0.022183547489474400, +0.007267287738381420,
            +0.000280675284928631, +0.004021769226522670, +0.004666323108506110, -0.007639556152399440, -0.023883316927296100,
            -0.027429491512297200, -0.015701400458900900, -0.002958012517850590, -0.000983630371410803, -0.004298625913755040,
            +0.000497235250104881, +0.014818029606027200, +0.025171758121561500, +0.020792218661758700, +0.007570543693995500,
            -0.000149653919325934, +0.001750736605729410, +0.003089076325810160, -0.005305308994708100, -0.018108455455086400,
            -0.022047333565259600, -0.013267551943819000, -0.002018363284070140, +0.001069639556348870, -0.001860209176511530,
            -0.000463614071435545, +0.008704543543861430, +0.017143085286423300, +0.015638408637062000, +0.006296408560066590,
            -0.000535925742056348, +0.000171657836788521, +0.002440508360430670, -0.001418894179803900, -0.009776930284869270,
            -0.013854882455851700, -0.009509493715276500, -0.002402759588650490, -0.000126769418915459, -0.002455878063427290
        };

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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(5e-11);
    }

    [TestMethod]
    public void TypeII_Even_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.014978759442673300, +0.024213566041449500, +0.022188255529808500, -0.011118960862351900, -0.052771796615212700,
            -0.080990391362511100, -0.059416835744030200, +0.018009279363915900, +0.107669861136537000, +0.144958782721866000,
            +0.094009385661312700, -0.021122818851714300, -0.128738485213722000, -0.159246506507874000, -0.096252501126483800,
            +0.012003241661359500, +0.091505253379345200, +0.099587968417279500, +0.052340878657204700, +0.001509352218966960,
            -0.014502171505733300, -0.000066915739404673, +0.011245496870296200, -0.004866090952279130, -0.039066233202330200,
            -0.059562900127856200, -0.044524429264181300, -0.002872839695333230, +0.035092115121871100, +0.045591979860382900,
            +0.030194823189362100, +0.009613926996574650, +0.001633331520038740, +0.005755542436792900, +0.007795207654586580,
            -0.002889288328895110, -0.021820368588420900, -0.034075411598412800, -0.029024608294319500, -0.010059478500377000,
            +0.009337592000632010, +0.018035992925089600, +0.015788021548616600, +0.010682339229426800, +0.009611627835045610,
            +0.011764306321616600, +0.010867695971714800, +0.002795543506369170, -0.009563602973777230, -0.018700185992591400,
            -0.019214859862975300, -0.012168906991686600, -0.003324962859813360, +0.002335603712577620, +0.004327316117023350,
            +0.005542960377646000, +0.008194999027986320, +0.011199936585246900, +0.011396447010625300, +0.006978848930884780,
            -0.000441605742923806, -0.007050968121008230, -0.009929744349876560, -0.009039938905440080, -0.006499090561129910,
            -0.004237239095878490, -0.002444494056096580, -9.84605455818806e-05, +0.003309791904375270, +0.006766991835553470,
            +0.008428940290460120, +0.007229624265207330, +0.003808240514366860, -5.45888262092609e-05, -0.002839383880416040,
            -0.004199726127331250, -0.004693247428208400, -0.004826643271776340, -0.004459058883598450, -0.003093618062839660,
            -0.000639021287652974, +0.002200299271199510, +0.004307562306678230, +0.004932093560538300, +0.004143227345328140,
            +0.002608371888095290, +0.000990977849630774, -0.000447910393354501, -0.001750804133498940, -0.002890320597565150,
            -0.003569851122642080, -0.003408708135062860, -0.002307790529783660, -0.000620292178357255, +0.001040233167470430,
            +0.002174606643479950, +0.002631124732550260, +0.002541032735334680, +0.002085216525012540, +0.001342755609884550
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod]
    public void TypeII_Odd_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.018631883268967100, +0.033242378741733600, +0.025793939719573800, -0.014547489430958800, -0.080881877450630700,
            -0.109077012595877000, -0.063683051969402300, +0.049501110342035800, +0.162026375489273000, +0.182326480269357000,
            +0.078879220643182900, -0.081909952641598600, -0.186256605004734000, -0.166140869133506000, -0.051553968352943200,
            +0.061993882347085800, +0.098313173846301800, +0.059870858088622200, +0.010291414015800200, +0.003314858664158320,
            +0.031686960128460900, +0.044649805263802800, +0.009512302725429230, -0.051064358540042100, -0.083945787411508100,
            -0.061860810584356800, -0.009315765694124060, +0.026696622639163200, +0.025452330103701100, +0.007823890519983960,
            +0.006137645801748450, +0.026165815297591900, +0.043118328277092200, +0.032317623393412800, -0.002385758583597520,
            -0.033245557694446500, -0.038331887224896800, -0.022210823334034400, -0.006786818929099280, -0.005678486468247660,
            -0.011477541197403700, -0.007865733819720940, +0.009689623369064270, +0.028384580968652200, +0.032341152421605100,
            +0.019343345055753400, +0.001831206077536050, -0.007215764802567310, -0.006821477678592540, -0.005993357715613420,
            -0.011344485500707200, -0.019183953313584500, -0.020272046438802700, -0.010704088953410000, +0.003610915862057740,
            +0.013182598752932500, +0.014214132873217800, +0.010857039389329100, +0.009051521629415150, +0.009689654101221560,
            +0.008561244882763530, +0.002261697895311400, -0.007175215814779680, -0.013993186999123500, -0.014421020126974700,
            -0.009936150853576890, -0.004942614199611080, -0.002008600741013500, +1.25620621392513e-05, +0.003440708469094610,
            +0.008241285381879060, +0.011605531454022900, +0.010868206169902300, +0.006302061712568230, +0.000782919153464336,
            -0.003032644384844520, -0.004819226572961160, -0.005942479099369400, -0.007216890403585610, -0.007768177760300830,
            -0.006148611010283620, -0.002230698010497550, +0.002320657304416670, +0.005451388077722360, +0.006393031208349630,
            +0.005887510459671300, +0.004969454929373870, +0.003820732405842630, +0.001908532820429470, -0.000949614092876955,
            -0.003939000233013840, -0.005769546869653470, -0.005761086233614030, -0.004322309109938610, -0.002403241617607810,
            -0.000614444663324253, +0.001050724975163450, +0.002736020379905210, +0.004171901667848610, +0.004729646654726940
        };

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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod]
    public void TypeII_Corrected_Even_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.027988053972296900, +0.062810840221972700, +0.057212037751087700, -0.032413524924819800, -0.155449712245133000,
            -0.201340089046208000, -0.081673207212271200, +0.122360788589317000, +0.239648442095990000, +0.178481468937220000,
            +0.016820632799034200, -0.086866172714557300, -0.068732968535150000, -0.010618063442382300, -0.013306747348527100,
            -0.068522705216921600, -0.090065978394327000, -0.039771381388358800, +0.028306148290808300, +0.048769564366496900,
            +0.026073260320450100, +0.012613418832042100, +0.031085012092532300, +0.051014261117942300, +0.039074708787122900,
            +0.004267364515620640, -0.019383164023757200, -0.018791067384440100, -0.012296484961069800, -0.018371916649107000,
            -0.030970880744940000, -0.031907600388687300, -0.017119594401438300, -0.000256477956192067, +0.006811438529769910,
            +0.007294472100851310, +0.010956994535936600, +0.019146544765477500, +0.023977910510257100, +0.020020911111180400,
            +0.011042133278308000, +0.004035087379887450, +0.000554063658137516, -0.003150167303346150, -0.009236415005702550,
            -0.014924187485088200, -0.016402496105707900, -0.013626288742930800, -0.009684048111074260, -0.006572895028375020,
            -0.003433075868223220, +0.001067315443497130, +0.006173371171857250, +0.009730794638366240, +0.010689767152999900,
            +0.009915497403778650, +0.008620546295284320, +0.006886818484718830, +0.004119693529588900, +0.000415584335887277,
            -0.003162436822515260, -0.005616400439738800, -0.006840478789520730, -0.007272927282160740, -0.007078783152482900,
            -0.006010566733276210, -0.003978075030179170, -0.001423568708404700, +0.000988182436638831, +0.002886175854489640,
            +0.004263785598642890, +0.005142645858307640, +0.005370161231978280, +0.004798329326979330, +0.003537462343457760,
            +0.001927166051533770, +0.000288137468662351, -0.001219483644736440, -0.002507330129167190, -0.003432288730843160,
            -0.003828381519410630, -0.003641368697461960, -0.002977827235270750, -0.002015654647824830, -0.000901641226357868,
            +0.000252788740826069, +0.001317387767456570, +0.002136996605183060, +0.002593872953726820, +0.002660705619384410,
            +0.002385877404102190, +0.001842679994235150, +0.001106763792112070, +0.000271443125875640, -0.000543631188453204,
            -0.001218837256794820, -0.001672492737539670, -0.001873069751020860, -0.001823706496663390, -0.001547222208480880
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.032046244824887300, +0.071627442221126000, +0.060799080936585900, -0.044935648045190500, -0.189632067746830000,
            -0.219497096142217000, -0.066695223546211800, +0.173081196083005000, +0.282179129693172000, +0.171936487965410000,
            -0.022446697385719100, -0.111747602592483000, -0.062110670666474000, -0.000531906728658865, -0.029716653391800800,
            -0.099724695816931000, -0.097090905058943000, -0.010744701295503700, +0.060196803545379100, +0.051199197847827000,
            +0.011518601591762000, +0.016242892467179400, +0.056606091606206000, +0.065228669313933400, +0.021292229174587800,
            -0.024313520835196600, -0.027927517982712900, -0.009972603193427500, -0.013577740972597900, -0.038647834651119500,
            -0.048113085845902900, -0.026260025824378800, +0.001279069010767870, +0.008225346998165540, +0.002504015613617340,
            +0.007483793771032050, +0.024825123672031500, +0.034354734895956400, +0.025361516861845200, +0.010293511911602200,
            +0.004324920942429720, +0.004834465232773820, -0.000537520990923501, -0.013391687773719900, -0.022669214745099000,
            -0.021215176365308400, -0.014664433497086400, -0.011425542952586400, -0.010849699907745400, -0.006421697695639710,
            +0.003062209213705850, +0.011492179457759300, +0.014124150719775600, +0.013162037902239100, +0.012991834722225900,
            +0.013483359754173400, +0.011089646967685900, +0.004897535715734340, -0.001762749497556640, -0.005814039168720590,
            -0.007854105040952720, -0.009977478402225940, -0.012050445384660200, -0.011984075844528000, -0.008999199505855320,
            -0.004792257483439110, -0.001246033112139870, +0.001618899431696250, +0.004700603099045530, +0.007778139663586220,
            +0.009500192209220670, +0.009148701781380410, +0.007442805356026290, +0.005406976783071380, +0.003178188087519690,
            +0.000389169577457936, -0.002728995299207660, -0.005261150013577360, -0.006562661679887330, -0.006798529804880230,
            -0.006410410020473250, -0.005445916864447090, -0.003705962980986550, -0.001342409398865180, +0.001050655635921650,
            +0.002951144100175120, +0.004261218778096470, +0.005074917926704170, +0.005316688006076580, +0.004799212495150650,
            +0.003554144453966710, +0.001912703378550080, +0.000235285032879509, -0.001310070773379400, -0.002657126594610970,
            -0.003654969650061180, -0.004095336801310540, -0.003902707202644900, -0.003202857885700300, -0.002189151222260010
        };

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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod]
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

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

        k_low_db.ToDebug();
        k_sl99_db.ToDebug();
        k_sl_db.ToDebug();

        k_sl_pl_db.ToDebug();
        k_pl_sl_db.ToDebug();

        k_pl_db.ToDebug();
        k_pl_ph_db.ToDebug();
        k_c0_db.ToDebug();
        k_ph_pl_db.ToDebug();
        k_ph_db.ToDebug();

        k_ph_sh_db.ToDebug();
        k_sh_ph_db.ToDebug();

        k_sh_db.ToDebug();
        k_sh_fd05_db.ToDebug();
        k_fd05_db.ToDebug();

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
        k_ph_db.AssertGreaterOrEqualsThan(-Rp, 6.13e-3);      // Коэффициент передачи на верхней границе полосы подавления

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

        h_pl.AssertGreaterOrEqualsThan(-Rp, 1.17e-4);
        h_pl_ph.AssertGreaterOrEqualsThan(-Rp);
        h_c0.AssertGreaterOrEqualsThan(-Rp);
        h_ph_pl.AssertGreaterOrEqualsThan(-Rp);
        h_ph.AssertGreaterOrEqualsThan(-Rp, 4.9e-6);

        h_ph_sh.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);

        h_sh.AssertLessOrEqualsThan(-Rs);
        h_sh_fd05.AssertLessOrEqualsThan(-Rs);
        h_fd05.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
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
        const double fph = 12.5 / Consts.pi2; // верхняя частота границы полосы заграждения
        const double fsh = 15 / Consts.pi2; // верхняя частота границы полосы пропускания

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

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

        k_low_db.ToDebug();
        k_sl99_db.ToDebug();
        k_sl_db.ToDebug();

        k_sl_pl_db.ToDebug();
        k_pl_sl_db.ToDebug();

        k_pl_db.ToDebug();
        k_pl_ph_db.ToDebug();
        k_c0_db.ToDebug();
        k_ph_pl_db.ToDebug();
        k_ph_db.ToDebug();

        k_ph_sh_db.ToDebug();
        k_sh_ph_db.ToDebug();

        k_sh_db.ToDebug();
        k_sh_fd05_db.ToDebug();
        k_fd05_db.ToDebug();

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
        k_ph_db.AssertGreaterOrEqualsThan(-Rp, 8.1e-3);      // Коэффициент передачи на верхней границе полосы подавления

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

        h_pl.AssertGreaterOrEqualsThan(-Rp, 4.56e-4);
        h_pl_ph.AssertGreaterOrEqualsThan(-Rp);
        h_c0.AssertGreaterOrEqualsThan(-Rp);
        h_ph_pl.AssertGreaterOrEqualsThan(-Rp);
        h_ph.AssertGreaterOrEqualsThan(-Rp, 4.9e-6);

        h_ph_sh.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);

        h_sh.AssertLessOrEqualsThan(-Rs);
        h_sh_fd05.AssertLessOrEqualsThan(-Rs);
        h_fd05.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

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
        var x_pl_ph = GetSinSignal(fpl + (fph - fpl) * 0.12);        // частота выше нижней границы пропускания
        var x_c0 = GetSinSignal((fpl * fph).Sqrt());                // частота в середине полосы пропускания
        var x_ph_pl = GetSinSignal(fph - (fph - fpl) * 0.25);        // частота ниже верхней границы пропускания
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

        k_low_db.ToDebug();
        k_sl99_db.ToDebug();
        k_sl_db.ToDebug();

        k_sl_pl_db.ToDebug();
        k_pl_sl_db.ToDebug();

        k_pl_db.ToDebug();
        k_pl_ph_db.ToDebug();
        k_c0_db.ToDebug();
        k_ph_pl_db.ToDebug();
        k_ph_db.ToDebug();

        k_ph_sh_db.ToDebug();
        k_sh_ph_db.ToDebug();

        k_sh_db.ToDebug();
        k_sh_fd05_db.ToDebug();
        k_fd05_db.ToDebug();

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
        k_pl_db.AssertLessOrEqualsThan(-Rs, 0.86);      // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_pl_ph_db.AssertGreaterOrEqualsThan(-Rp * 2);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertGreaterOrEqualsThan(-Rp, 1.5e-3);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_ph_pl_db.AssertGreaterOrEqualsThan(-Rp * 2);    // Коэффициент передачи у верхнего края полосы подавления
        k_ph_db.AssertLessOrEqualsThan(-Rs, 0.3);      // Коэффициент передачи на верхней границе полосы подавления

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
        var h_pl_ph = H.GetValue(fpl + (fph - fpl) * 0.12).Power.In_dB_byPower();
        var h_c0 = H.GetValue((fpl * fph).Sqrt()).Power.In_dB_byPower();
        var h_ph_pl = H.GetValue(fph - (fph - fpl) * 0.25).Power.In_dB_byPower();
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

        h_pl.AssertLessOrEqualsThan(-Rs);
        h_pl_ph.AssertGreaterOrEqualsThan(-Rp * 2);
        h_c0.AssertGreaterOrEqualsThan(-Rp, 1.56e-5);
        h_ph_pl.AssertGreaterOrEqualsThan(-Rp * 2);
        h_ph.AssertLessOrEqualsThan(-Rs);

        h_ph_sh.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);

        h_sh.AssertLessOrEqualsThan(-Rs);
        h_sh_fd05.AssertLessOrEqualsThan(-Rs);
        h_fd05.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
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
        const double fph = 12.5 / Consts.pi2; // верхняя частота границы полосы заграждения
        const double fsh = 15 / Consts.pi2; // верхняя частота границы полосы пропускания

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

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
        var x_pl_ph = GetSinSignal(fpl + (fph - fpl) * 0.12);        // частота выше нижней границы пропускания
        var x_c0 = GetSinSignal((fpl * fph).Sqrt());                // частота в середине полосы пропускания
        var x_ph_pl = GetSinSignal(fph - (fph - fpl) * 0.25);        // частота ниже верхней границы пропускания
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

        k_low_db.ToDebug();
        k_sl99_db.ToDebug();
        k_sl_db.ToDebug();

        k_sl_pl_db.ToDebug();
        k_pl_sl_db.ToDebug();

        k_pl_db.ToDebug();
        k_pl_ph_db.ToDebug();
        k_c0_db.ToDebug();
        k_ph_pl_db.ToDebug();
        k_ph_db.ToDebug();

        k_ph_sh_db.ToDebug();
        k_sh_ph_db.ToDebug();

        k_sh_db.ToDebug();
        k_sh_fd05_db.ToDebug();
        k_fd05_db.ToDebug();

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
        k_pl_db.AssertLessOrEqualsThan(-Rs, 2.36);      // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_pl_ph_db.AssertGreaterOrEqualsThan(-Rp * 2);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertGreaterOrEqualsThan(-Rp, 1.5e-3);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_ph_pl_db.AssertGreaterOrEqualsThan(-Rp * 2);    // Коэффициент передачи у верхнего края полосы подавления
        k_ph_db.AssertLessOrEqualsThan(-Rs, 1.57);      // Коэффициент передачи на верхней границе полосы подавления

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
        var h_pl_ph = H.GetValue(fpl + (fph - fpl) * 0.12).Power.In_dB_byPower();
        var h_c0 = H.GetValue((fpl * fph).Sqrt()).Power.In_dB_byPower();
        var h_ph_pl = H.GetValue(fph - (fph - fpl) * 0.25).Power.In_dB_byPower();
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

        h_pl.AssertLessOrEqualsThan(-Rs);
        h_pl_ph.AssertGreaterOrEqualsThan(-Rp * 2);
        h_c0.AssertGreaterOrEqualsThan(-Rp, 1.56e-5);
        h_ph_pl.AssertGreaterOrEqualsThan(-Rp * 2);
        h_ph.AssertLessOrEqualsThan(-Rs);

        h_ph_sh.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);

        h_sh.AssertLessOrEqualsThan(-Rs);
        h_sh_fd05.AssertLessOrEqualsThan(-Rs);
        h_fd05.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
    public void TypeII_Corrected_Even_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

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
        k_pl_ph_db.AssertGreaterOrEqualsThan(-Rp, 3.67e-3);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertGreaterOrEqualsThan(-Rp, 9.14e-4);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_ph_pl_db.AssertGreaterOrEqualsThan(-Rp, 4.06e-2);    // Коэффициент передачи у верхнего края полосы подавления
        k_ph_db.AssertGreaterOrEqualsThan(-Rp, 0.453);      // Коэффициент передачи на верхней границе полосы подавления

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

        h_pl.AssertGreaterOrEqualsThan(-Rp, 0.452);
        h_pl_ph.AssertGreaterOrEqualsThan(-Rp, 1.68e-3);
        h_c0.AssertGreaterOrEqualsThan(-Rp);
        h_ph_pl.AssertGreaterOrEqualsThan(-Rp, 3.92e-2);
        h_ph.AssertGreaterOrEqualsThan(-Rp, 0.452);

        h_ph_sh.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);

        h_sh.AssertLessOrEqualsThan(-Rs);
        h_sh_fd05.AssertLessOrEqualsThan(-Rs);
        h_fd05.AssertLessOrEqualsThan(-Rs);
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_SignalProcessing()
    {
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

        var filter = new DSP.Filters.ChebyshevBandPass(dt, fsl, fpl, fph, fsh, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

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
        k_pl_ph_db.AssertGreaterOrEqualsThan(-Rp, 3.67e-3);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertGreaterOrEqualsThan(-Rp, 9.14e-4);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_ph_pl_db.AssertGreaterOrEqualsThan(-Rp, 4.06e-2);    // Коэффициент передачи у верхнего края полосы подавления
        k_ph_db.AssertGreaterOrEqualsThan(-Rp, 0.453);      // Коэффициент передачи на верхней границе полосы подавления

        // Коэффициенты передачи в переходной полосе
        // между верхней полосой заграждения и верхней
        // полосой пропускания
        k_ph_sh_db.AssertLessThan(-Rp, 0.396);
        k_sh_ph_db.AssertLessOrEqualsThan(-Rp);
        k_sh_ph_db.AssertLessOrEqualsThan(-Rp);

        // Коэффициенты передачи в верхней полосе пропускания
        // от верхней границы полосы пропускания до частоты,
        // близкой к половине частоты дискретизации
        k_sh_db.AssertLessOrEqualsThan(-Rs, 1.16);
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

        h_pl.AssertGreaterOrEqualsThan(-Rp, 0.452);
        h_pl_ph.AssertGreaterOrEqualsThan(-Rp, 1.68e-3);
        h_c0.AssertGreaterOrEqualsThan(-Rp);
        h_ph_pl.AssertGreaterOrEqualsThan(-Rp, 3.92e-2);
        h_ph.AssertGreaterOrEqualsThan(-Rp, 0.452);

        h_ph_sh.AssertLessOrEqualsThan(-Rp, 0.399);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_sh_ph.AssertLessOrEqualsThan(-Rp);

        h_sh.AssertLessOrEqualsThan(-Rs, 7.37e-4);
        h_sh_fd05.AssertLessOrEqualsThan(-Rs);
        h_fd05.AssertLessOrEqualsThan(-Rs);
    }
}