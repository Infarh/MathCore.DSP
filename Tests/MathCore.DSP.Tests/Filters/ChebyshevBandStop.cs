using System.Diagnostics;
using System.Globalization;
using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;
using static MathCore.Polynom.Array;
using System.Linq;
using MathCore.DSP.Signals;
// ReSharper disable RedundantArgumentDefaultValue

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ChebyshevBandStop : ChebyshevFiltersTests
{
    [TestMethod]
    public void TypeI_Creation()
    {
        const double fd = 10;           // Частота дискретизации
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

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
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

        // Выбор опорной частоты
        // Если   Wc / Wph > Wpl
        // то есть      Wc > Wpl*Wph
        // то есть Wsl*Wsh > Wpl*Wph
        // то есть центральная частота по границам подавления > центральной частоты по границам пропускания
        // то выбираем в качестве опорной частоты выбираем верхнюю границу пропускания
        // иначе, выбираем нижнюю границу пропускания
        var Wp = Wc / Wph > Wpl
            ? Wph
            : Wpl;
        var W0 = Abs(dW * Wp / (Wc - Wp.Pow2()));   // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F0 = W0 / Consts.pi2;
        const double F1 = 1 / Consts.pi2;

        W0.AssertEquals(0.615059351152204);
        F0.AssertEquals(0.0978897360307671);

        var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.50884713990958752);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.51284645671973);
        kW.AssertEquals(1.6258593550796006);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(6);

        var beta = arcsh(1 / eps_p) / N;
        beta.AssertEquals(0.2379958931439209);

        var sh = Sinh(beta) * W0;
        var ch = Cosh(beta) * W0;

        var poles = new Complex[N];
        if (N.IsOdd()) poles[0] = -sh;
        var r = N % 2;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var th = dth * (i - r + 1);

            var sin = -sh * Sin(th);
            var cos = ch * Cos(th);
            //poles[i] = new Complex(sin, cos);
            //poles[i + 1] = poles[i].ComplexConjugate;
            (poles[i], poles[i + 1]) = Complex.Conjugate(sin, cos);
        }

        poles.AssertEquals(
            (-0.038245020148109342, 0.611006849626110071),
            (-0.038245020148109342, -0.611006849626110071),
            (-0.104487338181130124, 0.447288057698909503),
            (-0.104487338181130124, -0.447288057698909503),
            (-0.142732358329239445, 0.163718791927200735),
            (-0.142732358329239445, -0.163718791927200735));

        var sqrtWc = Wc.Sqrt();
        var pzf_zeros = Enumerable.Range(0, 2 * N).Select(i => i % 2 == 0 ? Complex.ImValue(sqrtWc) : Complex.ImValue(-sqrtWc));
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh).ToArray();

        pzf_poles.ToDebugEnum();

        pzf_zeros.AssertEquals(
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235));

        pzf_poles.AssertEquals(
            (-0.134740876884887117, 2.966182443413979186),
            (-0.847786715984386374, -18.663156503433079081),
            (-0.134740876884887117, -2.966182443413979186),
            (-0.847786715984386374, 18.663156503433079081),
            (-0.446533081604365112, 2.352017474088679450),
            (-4.321880813163225810, -22.764582541042127417),
            (-0.446533081604365112, -2.352017474088679450),
            (-4.321880813163225810, 22.764582541042127417),
            (-0.792057068460424674, 0.960760828695931934),
            (-28.339128814192328321, -34.375205989850492472),
            (-0.792057068460424674, -0.960760828695931934),
            (-28.339128814192328321, 34.375205989850492472));

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        //z_zeros.Foreach(DWL);

        z_zeros.AssertEquals(
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148));

        //z_poles.Foreach(DWL);
        z_poles.AssertEquals(
            (0.944417945958803351, 0.286445125329786832),
            (0.065097723009304609, -0.953486610678531687),
            (0.944417945958803351, -0.286445125329786832),
            (0.065097723009304609, 0.953486610678531687),
            (0.930772935510864707, 0.222101794211007481),
            (-0.123362161088371003, -0.820507862685686762),
            (0.930772935510864707, -0.222101794211007481),
            (-0.123362161088371003, 0.820507862685686762),
            (0.919712629269422322, 0.088706215574633693),
            (-0.450430385390868437, -0.390813181192564196),
            (0.919712629269422322, -0.088706215574633693),
            (-0.450430385390868437, 0.390813181192564196));

        var G_norm = (N.IsOdd() ? 1 : Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        G_norm.AssertEquals(0.0342439351432325);

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var h_f00 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0, dt);
        var h_fpl = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fpl, dt);
        var h_fsl = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fsl, dt);
        var h_fcc = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, (fsl * fsh).Sqrt(), dt);
        var h_fsh = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fsh, dt);
        var h_fph = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fph, dt);
        var h_fd5 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fd / 2, dt);

        //var h_f00_db = h_f00.Power.In_dB_byPower().ToDebug();
        //var h_fpl_db = h_fpl.Power.In_dB_byPower().ToDebug();
        //var h_fsl_db = h_fsl.Power.In_dB_byPower().ToDebug();
        //var h_fcc_db = h_fcc.Power.In_dB_byPower().ToDebug();
        //var h_fsh_db = h_fsh.Power.In_dB_byPower().ToDebug();
        //var h_fph_db = h_fph.Power.In_dB_byPower().ToDebug();
        //var h_fd5_db = h_fd5.Power.In_dB_byPower().ToDebug();

        h_f00.Abs.AssertGreaterOrEqualsThan(Gp, 6.3e-12);
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fcc.Abs.AssertLessOrEqualsThan(Gs);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp, 1.5e-14);
        h_fd5.Abs.AssertGreaterOrEqualsThan(Gp, 7.5e-15);

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

        Assert.That.Collection(filter.A).IsEqualTo(A, 7.2e-15);
        Assert.That.Collection(filter.B).IsEqualTo(B, 5.6e-14);
    }

    [TestMethod]
    public void TypeII_Creation()
    {
        const double fd = 10;           // Частота дискретизации
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

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
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

        // Выбор опорной частоты
        // Если   Wc / Wph > Wpl
        // то есть      Wc > Wpl*Wph
        // то есть Wsl*Wsh > Wpl*Wph
        // то есть центральная частота по границам подавления > центральной частоты по границам пропускания
        // то выбираем в качестве опорной частоты выбираем верхнюю границу пропускания
        // иначе, выбираем нижнюю границу пропускания
        var Wp = Wc / Wph > Wpl
            ? Wph
            : Wpl;
        var W0 = Abs(dW * Wp / (Wc - Wp.Pow2()));   // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F0 = W0 / Consts.pi2;
        const double F1 = 1 / Consts.pi2;

        W0.AssertEquals(0.615059351152204);
        F0.AssertEquals(0.0978897360307671);

        var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.50884713990958752);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.51284645671973);
        kW.AssertEquals(1.6258593550796006);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(6);

        var beta = arcsh(1 / eps_p) / N;
        beta.AssertEquals(0.2379958931439209);

        var sh = Sinh(beta) * W0;
        var ch = Cosh(beta) * W0;

        var poles = new Complex[N];
        var r = N % 2;
        var shb = Sinh(beta);
        var chb = Cosh(beta);
        if (N.IsOdd()) poles[0] = -1 / shb;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var th = dth * (i - r + 1);
            var sin = Sin(th);
            var cos = Cos(th);
            var norm = 1 / (sin * sin * shb * shb + cos * cos * chb * chb);
            (poles[i], poles[i + 1]) = Complex.Conjugate(-shb * sin * norm, chb * cos * norm);
        }

        var L = N / 2;
        var zeros = new Complex[2 * L];
        for (var (n, dth) = (1, PI / N); n <= L; n++)
        {
            var th = dth * (n - 0.5);
            (zeros[2 * n - 2], zeros[2 * n - 1]) = Complex.Conjugate(0, 1 / Cos(th));
        }

        poles.AssertEquals(
            (-0.038245020148109342, 0.611006849626110071),
            (-0.038245020148109342, -0.611006849626110071),
            (-0.104487338181130124, 0.447288057698909503),
            (-0.104487338181130124, -0.447288057698909503),
            (-0.142732358329239445, 0.163718791927200735),
            (-0.142732358329239445, -0.163718791927200735));

        var sqrtWc = Wc.Sqrt();
        var pzf_zeros = Enumerable.Range(0, 2 * N).Select(i => i % 2 == 0 ? Complex.ImValue(sqrtWc) : Complex.ImValue(-sqrtWc));
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh).ToArray();

        pzf_poles.ToDebugEnum();

        pzf_zeros.AssertEquals(
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235),
            (0, 7.447990244669614235),
            (0, -7.447990244669614235));

        pzf_poles.AssertEquals(
            (-0.134740876884887117, 2.966182443413979186),
            (-0.847786715984386374, -18.663156503433079081),
            (-0.134740876884887117, -2.966182443413979186),
            (-0.847786715984386374, 18.663156503433079081),
            (-0.446533081604365112, 2.352017474088679450),
            (-4.321880813163225810, -22.764582541042127417),
            (-0.446533081604365112, -2.352017474088679450),
            (-4.321880813163225810, 22.764582541042127417),
            (-0.792057068460424674, 0.960760828695931934),
            (-28.339128814192328321, -34.375205989850492472),
            (-0.792057068460424674, -0.960760828695931934),
            (-28.339128814192328321, 34.375205989850492472));

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        //z_zeros.Foreach(DWL);

        z_zeros.AssertEquals(
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148),
            (0.756417559622531210, 0.654089042481751148),
            (0.756417559622531210, -0.654089042481751148));

        //z_poles.Foreach(DWL);
        z_poles.AssertEquals(
            (0.944417945958803351, 0.286445125329786832),
            (0.065097723009304609, -0.953486610678531687),
            (0.944417945958803351, -0.286445125329786832),
            (0.065097723009304609, 0.953486610678531687),
            (0.930772935510864707, 0.222101794211007481),
            (-0.123362161088371003, -0.820507862685686762),
            (0.930772935510864707, -0.222101794211007481),
            (-0.123362161088371003, 0.820507862685686762),
            (0.919712629269422322, 0.088706215574633693),
            (-0.450430385390868437, -0.390813181192564196),
            (0.919712629269422322, -0.088706215574633693),
            (-0.450430385390868437, 0.390813181192564196));

        var G_norm = (N.IsOdd() ? 1 : Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        G_norm.AssertEquals(0.0342439351432325);

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var h_f00 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0, dt);
        var h_fpl = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fpl, dt);
        var h_fsl = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fsl, dt);
        var h_fcc = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, (fsl * fsh).Sqrt(), dt);
        var h_fsh = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fsh, dt);
        var h_fph = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fph, dt);
        var h_fd5 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fd / 2, dt);

        //var h_f00_db = h_f00.Power.In_dB_byPower().ToDebug();
        //var h_fpl_db = h_fpl.Power.In_dB_byPower().ToDebug();
        //var h_fsl_db = h_fsl.Power.In_dB_byPower().ToDebug();
        //var h_fcc_db = h_fcc.Power.In_dB_byPower().ToDebug();
        //var h_fsh_db = h_fsh.Power.In_dB_byPower().ToDebug();
        //var h_fph_db = h_fph.Power.In_dB_byPower().ToDebug();
        //var h_fd5_db = h_fd5.Power.In_dB_byPower().ToDebug();

        //h_f00.Abs.AssertGreaterOrEqualsThan(Gp, 6.3e-12);
        //h_fpl.Abs.AssertGreaterOrEqualsThan(Gp);
        //h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        //h_fcc.Abs.AssertLessOrEqualsThan(Gs);
        //h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        //h_fph.Abs.AssertGreaterOrEqualsThan(Gp, 1.5e-14);
        //h_fd5.Abs.AssertGreaterOrEqualsThan(Gp, 7.5e-15);

        //var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        //Assert.That.Collection(filter.A).IsEqualTo(A, 7.2e-15);
        //Assert.That.Collection(filter.B).IsEqualTo(B, 5.6e-14);
    }

    [TestMethod]
    public void TypeI_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.03424393521213295, -0.15425500217716298, +0.35037692305779766, -0.41058062203391765, +0.20038885060936049,
            +0.14267212638406190, -0.18329409834716540, -0.00137539884933197, +0.14807095316906946, +0.04632156203643948,
            -0.08108247673066242, +0.03438061423219946, +0.15006589017607772, +0.05283773359118946, -0.01654832775770104,
            +0.07256759936603888, +0.13744615202055135, +0.08114563446737655, +0.03263544064249174, +0.07407374376606911,
            +0.11676644845322476, +0.08675317602309777, +0.03860846023149406, +0.04208768529398528, +0.06994236221495467,
            +0.05543360370915395, +0.00661472135001709, -0.01130905815867917, +0.00939906662610650, +0.00839277185320588,
            -0.03132518491630235, -0.05313587784868443, -0.03188307266991324, -0.01771455730697077, -0.04029162286261279,
            -0.05598314852023870, -0.03238562916654220, -0.00679143534523551, -0.01406394990151653, -0.02491948487338355,
            -0.00514501671814100, +0.02152842559479147, +0.02011670887801135, +0.00700210580042768, +0.01524586599216922,
            +0.03364589077558602, +0.03051186917817278, +0.01193748549310340, +0.00778587392704859, +0.01692309843292893,
            +0.01268134968903343, -0.00663859124216498, -0.01609205861754172, -0.00998532794760788, -0.00957185497879725,
            -0.02221649487437484, -0.02900036155812918, -0.02044054226260551, -0.01253733167580534, -0.01621365548677292,
            -0.01861940450416175, -0.00859815158203760, +0.00294641078794772, +0.00373333517174973, +0.00146378205416277,
            +0.00774699870278957, +0.01654277040711972, +0.01608490333855738, +0.01014629834640561, +0.00991905297325926,
            +0.01353296696374848, +0.01079248302453481, +0.00246838060989925, -0.00155728968537210, -0.00019328543210079,
            -0.00181877583802404, -0.00806455978523677, -0.01111698791150806, -0.00820680240464437, -0.00597537760805732,
            -0.00781595449010245, -0.00823103683899355, -0.00364272277241323, +0.00107101489791764, +0.00163063806415321,
            +0.00153937203074641, +0.00478397801822251, +0.00835462639434007, +0.00796491315123512, +0.00564828351736021,
            +0.00553282196908543, +0.00635338441554788, +0.00432738044610140, +0.00035161941407748, -0.00171081811405903,
            -0.00185475958466330, -0.00328607646040225, -0.00613299233598121, -0.00724343017879006, -0.00592237211590504
        };

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();

        Assert.That.Value(error2).LessThan(6e-15);
    }

    [TestMethod]
    public void TypeI_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        //const double total_time = 1 / fpl;
        //const int samples_count = (int)(total_time * fd) + 1;

        // Метод формирования гармонического сигнала на заданной частоте
        static EnumerableSignal GetSinSignal(double f0) => MathEnumerableSignal.Sin(dt, f0, (int)(100 * fd / (fpl * 0.1)));

        var x_low = GetSinSignal(fpl * 0.1);                    // частота, близкая к нулю
        var x_pl99 = GetSinSignal(fpl * 0.99);                   // частота чуть ниже нижней границы пропускания
        var x_pl = GetSinSignal(fpl);                          // частота на нижней границе пропускания

        var x_pl_sl = GetSinSignal(fpl + (fsl - fpl) * 0.1);      // частота выше нижней границы пропускания
        var x_sl_pl = GetSinSignal(fsl - (fsl - fpl) * 0.1);      // частота ниже нижней границы подавления

        var x_sl = GetSinSignal(fsl);                          // частота на границе подавления
        var x_sl_sh = GetSinSignal(fsl + (fsh - fsl) * 0.1);      // частота выше нижней границы подавления
        var x_c0 = GetSinSignal((fsl * fsh).Sqrt());           // частота в середине полосы подавления
        var x_sh_sl = GetSinSignal(fsh - (fsh - fsl) * 0.1);      // частота ниже верхней границы подавления
        var x_sh = GetSinSignal(fsh);                          // частота на границе подавления

        var x_sh_ph = GetSinSignal(fsh + (fph - fsh) * 0.1);      // частота выше верхней границы подавления
        var x_ph_sh = GetSinSignal(fph - (fph - fsh) * 0.1);      // частота ниже верхней границы пропускания

        var x_ph = GetSinSignal(fph);                          // частота на верхней границе пропускания
        var x_ph_fd05 = GetSinSignal(fph + (fd / 2 - fph) * 0.1);   // частота выше верхней границы пропускания
        var x_fd05 = GetSinSignal(0.9 * (fd / 2));               // частота ниже половины частоты дискретизации

        // Индивидуальная фильтрация каждой частотной составляющей
        var y_low = filter.ProcessIndividual(x_low);
        var y_pl99 = filter.ProcessIndividual(x_pl99);
        var y_pl = filter.ProcessIndividual(x_pl);

        var y_pl_sl = filter.ProcessIndividual(x_pl_sl);
        var y_sl_pl = filter.ProcessIndividual(x_sl_pl);

        var y_sl = filter.ProcessIndividual(x_sl);
        var y_sl_sh = filter.ProcessIndividual(x_sl_sh);
        var y_c0 = filter.ProcessIndividual(x_c0);
        var y_sh_sl = filter.ProcessIndividual(x_sh_sl);
        var y_sh = filter.ProcessIndividual(x_sh);

        var y_sh_ph = filter.ProcessIndividual(x_sh_ph);
        var y_ph_sh = filter.ProcessIndividual(x_ph_sh);

        var y_ph = filter.ProcessIndividual(x_ph);
        var y_ph_fd05 = filter.ProcessIndividual(x_ph_fd05);
        var y_fd05 = filter.ProcessIndividual(x_fd05);

        // Рассчёт коэффициентов передачи по мощности
        var k_low = y_low.Power / x_low.Power;
        var k_pl99 = y_pl99.Power / x_pl99.Power;
        var k_pl = y_pl.Power / x_pl.Power;

        var k_pl_sl = y_pl_sl.Power / x_pl_sl.Power;
        var k_sl_pl = y_sl_pl.Power / x_sl_pl.Power;

        var k_sl = y_sl.Power / x_sl.Power;
        var k_sl_sh = y_sl_sh.Power / x_sl_sh.Power;
        var k_c0 = y_c0.Power / x_c0.Power;
        var k_sh_sl = y_sh_sl.Power / x_sh_sl.Power;
        var k_sh = y_sh.Power / x_sh.Power;

        var k_sh_ph = y_sh_ph.Power / x_sh_ph.Power;
        var k_ph_sh = y_ph_sh.Power / x_ph_sh.Power;

        var k_ph = y_ph.Power / x_ph.Power;
        var k_ph_fd05 = y_ph_fd05.Power / x_ph_fd05.Power;
        var k_fd05 = y_fd05.Power / x_fd05.Power;

        // Рассчёт коэффициентов передачи по мощности в логарифмическим масштабе
        var k_low_db = k_low.In_dB_byPower();
        var k_pl99_db = k_pl99.In_dB_byPower();
        var k_pl_db = k_pl.In_dB_byPower();

        var k_pl_sl_db = k_pl_sl.In_dB_byPower();
        var k_sl_pl_db = k_sl_pl.In_dB_byPower();

        var k_sl_db = k_sl.In_dB_byPower();
        var k_sl_sh_db = k_sl_sh.In_dB_byPower();
        var k_c0_db = k_c0.In_dB_byPower();
        var k_sh_sl_db = k_sh_sl.In_dB_byPower();
        var k_sh_db = k_sh.In_dB_byPower();

        var k_sh_ph_db = k_sh_ph.In_dB_byPower();
        var k_ph_sh_db = k_ph_sh.In_dB_byPower();

        var k_ph_db = k_ph.In_dB_byPower();
        var k_ph_fd05_db = k_ph_fd05.In_dB_byPower();
        var k_fd05_db = k_fd05.In_dB_byPower();

        // Сравнение коэффициентов передачи с заданными параметрами фильтрации
        k_low_db.AssertGreaterOrEqualsThan(-Rp);        // Коэффициенты передачи в нижней полосе пропускания
        k_pl99_db.AssertGreaterOrEqualsThan(-Rp);       // должны быть меньше, чем заданный уровень
        k_pl_db.AssertGreaterOrEqualsThan(-Rp);         // неравномерности АЧХ (допуск) Rp - не более -1дБ

        // Коэффициенты передачи в переходной полосе
        // между полосой пропускания и полосой заграждения
        k_pl_sl_db.AssertGreaterOrEqualsThan(-Rs);      // Коэффициент передачи у нижнего края переходной полосы не должен быть ниже уровня подавления Rs (-40дБ), но может быть всё ещё порядка уровня пропускания Rp (-1дБ)
        k_sl_pl_db.AssertLessOrEqualsThan(-Rp);         // Коэффициент передачи у верхнего края переходной полосы должен быть гарантировано меньше коэффициента пропускания Rp (-1дБ) и должен приближаться к уровню Rs (-40дБ)

        // Коэффициенты передачи в полосе заграждения
        // должны бытьниже уровня подавления Rs
        k_sl_db.AssertLessOrEqualsThan(-Rs, 2.16);      // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_sl_sh_db.AssertLessOrEqualsThan(-Rs);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertLessOrEqualsThan(-Rs, 1.5);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_sh_sl_db.AssertLessOrEqualsThan(-Rs, 1.5);    // Коэффициент передачи у верхнего края полосы подавления
        k_sh_db.AssertLessOrEqualsThan(-Rs, 0.36);      // Коэффициент передачи на верхней границе полосы подавления

        // Коэффициенты передачи в переходной полосе
        // между верхней полосой заграждения и верхней
        // полосой пропускания
        k_sh_ph_db.AssertLessThan(-Rp);
        k_ph_sh_db.AssertGreaterOrEqualsThan(-Rs);
        k_ph_sh_db.AssertLessOrEqualsThan(-Rp);

        // Коэффициенты передачи в верхней полосе пропускания
        // от верхней границы полосы пропускания до частоты,
        // близкой к половине частоты дискретизации
        k_ph_db.AssertGreaterOrEqualsThan(-Rp, 0.01);
        k_ph_fd05_db.AssertGreaterOrEqualsThan(-Rp);
        k_fd05_db.AssertGreaterOrEqualsThan(-Rp);

        // Суммарный сигнал
        var x =
            x_low +
            x_pl99 +
            x_pl +
            x_pl_sl +
            x_sl_pl +
            x_sl +
            x_sl_sh +
            x_c0 +
            x_sh_sl +
            x_sh +
            x_sh_ph +
            x_ph_sh +
            x_ph +
            x_ph_fd05 +
            x_fd05;

        // Фильтрация суммарного гармонического сигнала
        var y = filter.ProcessIndividual(x);

        var X = x.GetSpectrum();    // Спектр входного сигнала фильтра
        var Y = y.GetSpectrum();    // Спектр выходного сигнала фильтра

        var H = Y / X;              // Коэффициент передачи фильтра как отношение спектров на выходе и на входе

        // Извлекаем из спетра коэффициенты передачи на частотах гармонических составляющих исходного сигнала
        var h_low = H.GetValue(fpl * 0.1).Power.In_dB_byPower();
        var h_pl99 = H.GetValue(fpl * 0.99).Power.In_dB_byPower();
        var h_pl = H.GetValue(fpl).Power.In_dB_byPower();

        var h_pl_sl = H.GetValue(fpl + (fsl - fpl) * 0.1).Power.In_dB_byPower();
        var h_sl_pl = H.GetValue(fsl - (fsl - fpl) * 0.1).Power.In_dB_byPower();

        var h_sl = H.GetValue(fsl).Power.In_dB_byPower();
        var h_sl_sh = H.GetValue(fsl + (fsh - fsl) * 0.1).Power.In_dB_byPower();
        var h_c0 = H.GetValue((fsl * fsh).Sqrt()).Power.In_dB_byPower();
        var h_sh_sl = H.GetValue(fsh - (fsh - fsl) * 0.1).Power.In_dB_byPower();
        var h_sh = H.GetValue(fsh).Power.In_dB_byPower();

        var h_sh_ph = H.GetValue(fsh + (fph - fsh) * 0.1).Power.In_dB_byPower();
        var h_ph_sh = H.GetValue(fph - (fph - fsh) * 0.1).Power.In_dB_byPower();

        var h_ph = H.GetValue(fph).Power.In_dB_byPower();
        var h_ph_fd05 = H.GetValue(fph + (fd / 2 - fph) * 0.1).Power.In_dB_byPower();
        var h_fd05 = H.GetValue(0.9 * (fd / 2)).Power.In_dB_byPower();

        // Тест фактических коэффициентов передачи
        h_low.AssertGreaterOrEqualsThan(-Rp);
        h_pl99.AssertGreaterOrEqualsThan(-Rp);
        h_pl.AssertGreaterOrEqualsThan(-Rp);

        h_pl_sl.AssertGreaterOrEqualsThan(-Rs);
        h_sl_pl.AssertLessOrEqualsThan(-Rp).GreaterOrEqualsThan(-Rs);

        h_sl.AssertLessOrEqualsThan(-Rs);
        h_sl_sh.AssertLessOrEqualsThan(-Rs);
        h_c0.AssertLessOrEqualsThan(-Rs);
        h_sh_sl.AssertLessOrEqualsThan(-Rs);
        h_sh.AssertLessOrEqualsThan(-Rs);

        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_ph_sh.AssertLessOrEqualsThan(-Rp).GreaterOrEqualsThan(-Rs);

        h_ph.AssertGreaterOrEqualsThan(-Rp);
        h_ph_fd05.AssertGreaterOrEqualsThan(-Rp);
        h_fd05.AssertGreaterOrEqualsThan(-Rp);
    }
}
