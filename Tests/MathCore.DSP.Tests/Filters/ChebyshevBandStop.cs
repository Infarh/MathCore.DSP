using System.Diagnostics;
using System.Globalization;
using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;
using static MathCore.Polynom.Array;
using System.Linq;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ChebyshevBandStop : UnitTest
{
    private static double arcsh(double x) => Log(x + Sqrt(x * x + 1));
    private static double arcch(double x) => Log(x + Sqrt(x * x - 1));

    [TestMethod]
    public void TypeI_Creation()
    {
        static void DWL(Complex z)
        {
            FormattableString formattable_string = $"({z.Re:F18}, {z.Im:F18}),";
            Debug.WriteLine(formattable_string.ToString(CultureInfo.InvariantCulture));
        }

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        Assert.That.Collection(filter.A).IsEqualTo(A, 7.2e-15);
        Assert.That.Collection(filter.B).IsEqualTo(B, 5.6e-14);
    }
}
