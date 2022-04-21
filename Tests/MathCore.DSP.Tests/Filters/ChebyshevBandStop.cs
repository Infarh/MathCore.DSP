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
using MathCore.DSP.Tests.Infrastructure;
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
        if (N.IsOdd())
            poles[0] = -sh;
        var r = N % 2;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (cos, sin) = Complex.Exp(dth * (i - r + 1));
            (poles[i], poles[i + 1]) = Complex.Conjugate(-sh * sin, ch * cos);
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

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(0.8830487276017475);

        //var sh = Sinh(beta) * W0;
        //var ch = Cosh(beta) * W0;

        var poles = new Complex[N];
        var r = N % 2;
        var sh = Sinh(beta);
        var ch = Cosh(beta);
        if (N.IsOdd()) poles[0] = -1 / sh;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (cos, sin) = Complex.Exp(dth * (i - r + 1));
            var z = new Complex(-sh * sin, ch * cos);
            var norm = W0 / z.Power;
            (poles[i], poles[i + 1]) = Complex.Conjugate(z.Re * norm, z.Im * norm);
        }

        var L = N / 2;
        var zeros = new Complex[2 * L];
        for (var (n, dth) = (1, PI / N); n <= L; n++)
            (zeros[2 * n - 2], zeros[2 * n - 1]) = Complex.Conjugate(0, W0 / Cos(dth * (n - 0.5)));

        poles.AssertEquals(
            /*[ 0]*/ (-0.082345846788124857, +0.434100948566762623),
            /*[ 1]*/ (-0.082345846788124857, -0.434100948566762623),
            /*[ 2]*/ (-0.289712304973535617, +0.409230909475116644),
            /*[ 3]*/ (-0.289712304973535617, -0.409230909475116644),
            /*[ 4]*/ (-0.555651050003078129, +0.210308138452255156),
            /*[ 5]*/ (-0.555651050003078129, -0.210308138452255156)
        );

        zeros.AssertEquals(
            Complex.ImValue(+0.6367562957863577),
            Complex.ImValue(-0.6367562957863577),
            Complex.ImValue(+0.8698252760638427),
            Complex.ImValue(-0.8698252760638427),
            Complex.ImValue(+2.3764068479140414),
            Complex.ImValue(-2.3764068479140414)
        );

        var pzf_zeros = AnalogBasedFilter.TransformToBandStop(zeros, Fsl, Fsh).ToArray();
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh).ToArray();

        pzf_zeros.AssertEquals(
            /*[ 0]*/ (+6.498578255002484E-16, +03.052371006978922274),
            /*[ 1]*/ (-6.498578255002484E-16, -18.173596380604330136),
            /*[ 2]*/ (+6.498578255002484E-16, +18.173596380604330136),
            /*[ 3]*/ (-6.498578255002484E-16, -03.052371006978922274),
            /*[ 4]*/ (+5.681953532024858E-16, +03.744581689426437876),
            /*[ 5]*/ (-5.681953532024858E-16, -14.814086935619917185),
            /*[ 6]*/ (+5.681953532024858E-16, +14.814086935619917185),
            /*[ 7]*/ (-5.681953532024858E-16, -03.744581689426437876),
            /*[ 8]*/ (+4.726274700415578E-16, +05.692732274427887873),
            /*[ 9]*/ (-4.726274700415578E-16, -09.744452401859820867),
            /*[10]*/ (+4.726274700415578E-16, +09.744452401859820867),
            /*[11]*/ (-4.726274700415578E-16, -05.692732274427887873)
        );

        pzf_poles.AssertEquals(
            /*[ 0]*/ (-00.357233675636420722, +02.285240672570521880),
            /*[ 1]*/ (-03.704090692484245828, -23.695242869460138024),
            /*[ 2]*/ (-00.357233675636420722, -02.285240672570521880),
            /*[ 3]*/ (-03.704090692484245828, +23.695242869460138024),
            /*[ 4]*/ (-01.289549889250637627, +02.373162913081646508),
            /*[ 5]*/ (-09.806194630472592877, -18.046372311366127406),
            /*[ 6]*/ (-01.289549889250637627, -02.373162913081646508),
            /*[ 7]*/ (-09.806194630472592877, +18.046372311366127406),
            /*[ 8]*/ (-03.239200436367614344, +02.141183612997055619),
            /*[ 9]*/ (-11.917868451476278935, -07.877976411602194418),
            /*[10]*/ (-03.239200436367614344, -02.141183612997055619),
            /*[11]*/ (-11.917868451476278935, +07.877976411602194418)
        );

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (0.954475531310739145, +0.298289222281131250),
            /*[ 1]*/ (0.095470822516272452, -0.995432228756968795),
            /*[ 2]*/ (0.095470822516272452, +0.995432228756968795),
            /*[ 3]*/ (0.954475531310739145, -0.298289222281131250),
            /*[ 4]*/ (0.932264972945861747, +0.361776201840657230),
            /*[ 5]*/ (0.291453286040585902, -0.956585062634862759),
            /*[ 6]*/ (0.291453286040585902, +0.956585062634862759),
            /*[ 7]*/ (0.932264972945861747, -0.361776201840657230),
            /*[ 8]*/ (0.850107950597354556, +0.526608462077059802),
            /*[ 9]*/ (0.616310667393904277, -0.787503118251908996),
            /*[10]*/ (0.616310667393904277, +0.787503118251908996),
            /*[11]*/ (0.850107950597354556, -0.526608462077059802)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (+0.940450693226825618, +0.217829048776244150),
            /*[ 1]*/ (-0.155948781938165276, -0.843736166288806211),
            /*[ 2]*/ (+0.940450693226825618, -0.217829048776244150),
            /*[ 3]*/ (-0.155948781938165276, +0.843736166288806211),
            /*[ 4]*/ (+0.855796430467403391, +0.206867091409867393),
            /*[ 5]*/ (-0.017983158020384617, -0.594569074184266189),
            /*[ 6]*/ (+0.855796430467403391, -0.206867091409867393),
            /*[ 7]*/ (-0.017983158020384617, +0.594569074184266189),
            /*[ 8]*/ (+0.706740793793766198, +0.157253491974098408),
            /*[ 9]*/ (+0.181254230329152499, -0.291557469659536383),
            /*[10]*/ (+0.706740793793766198, -0.157253491974098408),
            /*[11]*/ (+0.181254230329152499, +0.291557469659536383)
        );

        var zeros_to_poles = z_zeros.Multiply(z => 1 - z).Re / z_poles.Multiply(z => 1 - z).Re;

        var G_norm = 1 / zeros_to_poles;

        G_norm.AssertEquals(0.10613121177849884);

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        //A.ToDebugEnum();

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

        h_f00.Abs.AssertGreaterOrEqualsThan(Gp, 1e-11);
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp, 4.1e-2);

        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fcc.Abs.AssertLessOrEqualsThan(Gs);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);

        h_fph.Abs.AssertLessOrEqualsThan(Gp);
        h_fd5.Abs.AssertGreaterOrEqualsThan(Gp, 8.89e-16);

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        //B.ToDebugEnum();
        //A.ToDebugEnum();

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);
    }

    [TestMethod]
    public void TypeIICorrected_Creation()
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
        var kW = F0 / F1;

        kEps.AssertEquals(196.51284645671973);
        kW.AssertEquals(0.615059351152204);

        var N = (int)Ceiling(arcch(kEps) / arcch(1 / kW)); // Порядок фильтра
        N.AssertEquals(6);
        //N.ToDebug();

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(0.8830487276017475);
        //beta.ToDebug();

        //var sh = Sinh(beta) * W0;
        //var ch = Cosh(beta) * W0;

        var r = N % 2;
        var sh = Sinh(beta);
        var ch = Cosh(beta);

        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -(W0 / kW) / sh;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (cos, sin) = Complex.Exp(dth * (i - r + 1));
            var z = new Complex(-sin * sh, cos * ch);
            var norm = (W0 / kW) / z.Power;
            (poles[i], poles[i + 1]) = Complex.Conjugate(z.Re * norm, z.Im * norm);
        }
        //(W0 / kW).ToDebug();
        //poles.ToDebugEnum();

        var zeros = new Complex[N - r];
        for (var (n, dth, L) = (1, PI / N, N / 2); n <= L; n++)
        {
            var th = dth * (n - 0.5);
            (zeros[2 * n - 2], zeros[2 * n - 1]) = Complex.Conjugate(0, W0 / Cos(th) / kW);
        }

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.133882765352424299, 0.705787088276199515),
            /*[ 1]*/ (-0.133882765352424299, -0.705787088276199515),
            /*[ 2]*/ (-0.471031461322897249, 0.665351902557851593),
            /*[ 3]*/ (-0.471031461322897249, -0.665351902557851593),
            /*[ 4]*/ (-0.903410457807307554, 0.341931454351974973),
            /*[ 5]*/ (-0.903410457807307554, -0.341931454351974973)
        );

        //zeros.ToDebugEnum();
        zeros.AssertEquals(
            /*[ 0]*/ (0, +1.035276180410082958),
            /*[ 1]*/ (0, -1.035276180410082958),
            /*[ 2]*/ (0, +1.414213562373094923),
            /*[ 3]*/ (0, -1.414213562373094923),
            /*[ 4]*/ (0, +3.863703305156270140),
            /*[ 5]*/ (0, -3.863703305156270140)
        );

        var pzf_zeros = AnalogBasedFilter.TransformToBandStop(zeros, Fsl, Fsh).ToArray();
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh).ToArray();

        //pzf_zeros.ToDebugEnum();
        pzf_zeros.AssertEquals(AccuracyComplex.Eps(1e-15),
            /*[ 0]*/ (0, +04.130273096840153535),
            /*[ 1]*/ (0, -13.430724163768442025),
            /*[ 2]*/ (0, +13.430724163768442025),
            /*[ 3]*/ (0, -04.130273096840153535),
            /*[ 4]*/ (0, +04.784885622771338554),
            /*[ 5]*/ (0, -11.593288337071019356),
            /*[ 6]*/ (0, +11.593288337071019356),
            /*[ 7]*/ (0, -04.784885622771338554),
            /*[ 8]*/ (0, +06.305474677069272893),
            /*[ 9]*/ (0, -08.797523029697883246),
            /*[10]*/ (0, +08.797523029697883246),
            /*[11]*/ (0, -06.305474677069272893)
        );

        //pzf_poles.ToDebugEnum();
        pzf_poles.AssertEquals(
            /*[ 0]*/ (-0.418083970685262929, 3.313002346197622572),
            /*[ 1]*/ (-2.079871559989669549, -16.481424405583815940),
            /*[ 2]*/ (-0.418083970685262929, -3.313002346197622572),
            /*[ 3]*/ (-2.079871559989669549, 16.481424405583815940),
            /*[ 4]*/ (-1.458679273005567101, 3.598910399875938637),
            /*[ 5]*/ (-5.365862151846025441, -13.238864402857418057),
            /*[ 6]*/ (-1.458679273005567101, -3.598910399875938637),
            /*[ 7]*/ (-5.365862151846025441, 13.238864402857418057),
            /*[ 8]*/ (-3.337785088918605148, 4.449412112440246148),
            /*[ 9]*/ (-5.984711866607919539, -7.977880168844960096),
            /*[10]*/ (-3.337785088918605148, -4.449412112440246148),
            /*[11]*/ (-5.984711866607919539, 7.977880168844960096)
        );

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (0.918193111864285538, 0.396133070223857686),
            /*[ 1]*/ (0.378396915501328979, -0.925643438014379316),
            /*[ 2]*/ (0.378396915501328979, 0.925643438014379316),
            /*[ 3]*/ (0.918193111864285538, -0.396133070223857686),
            /*[ 4]*/ (0.891721948180851154, 0.452583657606577094),
            /*[ 5]*/ (0.496993846052344812, -0.867754064805286940),
            /*[ 6]*/ (0.496993846052344812, 0.867754064805286940),
            /*[ 7]*/ (0.891721948180851154, -0.452583657606577094),
            /*[ 8]*/ (0.819178186609497039, 0.573539099437149202),
            /*[ 9]*/ (0.675756207627617456, -0.737125191438157579),
            /*[10]*/ (0.675756207627617456, 0.737125191438157579),
            /*[11]*/ (0.819178186609497039, -0.573539099437149202)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.908793467740180883, 0.309717466443433231),
            /*[ 1]*/ (0.163387293585929599, -0.868405401795834253),
            /*[ 2]*/ (0.908793467740180883, -0.309717466443433231),
            /*[ 3]*/ (0.163387293585929599, 0.868405401795834253),
            /*[ 4]*/ (0.813050459740102038, 0.304073054638865137),
            /*[ 5]*/ (0.239332246978127994, -0.646828066383620781),
            /*[ 6]*/ (0.813050459740102038, -0.304073054638865137),
            /*[ 7]*/ (0.239332246978127994, 0.646828066383620781),
            /*[ 8]*/ (0.653844028518078435, 0.315309855864144728),
            /*[ 9]*/ (0.406761738442774023, -0.431906908686362612),
            /*[10]*/ (0.653844028518078435, -0.315309855864144728),
            /*[11]*/ (0.406761738442774023, 0.431906908686362612)
        );

        var re_to_im = z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z);

        var G_norm = 1 / re_to_im.Abs;

        G_norm.AssertEquals(0.21872794262004155);

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

        h_f00.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp);

        h_fsl.Abs.AssertLessOrEqualsThan(Gs, 1e-10);
        h_fcc.Abs.AssertLessOrEqualsThan(Gs);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs, 1e-12);

        h_fph.Abs.AssertGreaterOrEqualsThan(Gp, 1.5e-14);
        h_fd5.Abs.AssertGreaterOrEqualsThan(Gp, 7.5e-15);

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

        filter.A.AssertEquals(A);
        filter.B.AssertEquals(B);
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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();
        error2.ToDebug();

        Assert.That.Value(error2).LessThan(2e-19);
    }

    [TestMethod]
    public void TypeII_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.10613121177841096, -0.26103460227554337, +0.42236712570109500, -0.31236235334897255, +0.09459041763295879,
            +0.18016239552482743, -0.09916586355256542, +0.00930758408098219, +0.15853725937965005, +0.02593876266197550,
            -0.01409342485430792, +0.11046738664186878, +0.10078325785164347, +0.01659795072560020, +0.05870357651937367,
            +0.11251312120622058, +0.06538763967889322, +0.03895494438512425, +0.07798559977247310, +0.07820001634527280,
            +0.03948731841364082, +0.03930900904330541, +0.05395684424027589, +0.03413188140621725, +0.01314213394563154,
            +0.01740017749052419, +0.01461064011612873, -0.00414461262894429, -0.01185079618584469, -0.00993415959983337,
            -0.01744822109535300, -0.02706335674372616, -0.02676761425524501, -0.02551505689351452, -0.02926149191211634,
            -0.02956249001558693, -0.02477414466158015, -0.02183171830160506, -0.02007789199184376, -0.01498329574144289,
            -0.00887383178691553, -0.00482463814723077, -0.00064401562219056, +0.00477825177168189, +0.00904427331792642,
            +0.01169883003958944, +0.01433270456309558, +0.01655027770634429, +0.01715917082069593, +0.01677092419963938,
            +0.01611210341768801, +0.01462986475691491, +0.01215965079493991, +0.00941773925664902, +0.00660356277236413,
            +0.00345535138840489, +0.00024371567513888, -0.00259666710558731, -0.00509963460425009, -0.00731281727536239,
            -0.00898466278664610, -0.01000450782925266, -0.01049998839050708, -0.01048347936628878, -0.00988818853468350,
            -0.00881178661467246, -0.00740596391144897, -0.00572172074681069, -0.00382785492895453, -0.00187620130423637,
            +0.00001763467584617, +0.00179617670603852, +0.00337536278992109, +0.00466318259323689, +0.00562266266162753,
            +0.00624281850053863, +0.00650346813195885, +0.00640755862538227, +0.00599494265983126, +0.00531046403542551,
            +0.00439849207043666, +0.00332175423706104, +0.00215242692078924, +0.00095370668772464, -0.00021511407305229,
            -0.00129342019903346, -0.00223061053535510, -0.00299134013682365, -0.00354924733241965, -0.00388743061474759,
            -0.00400393748145245, -0.00390940040144447, -0.00362221666521721, -0.00316921995594049, -0.00258559454056941,
            -0.00191036263714600, -0.00118343873029951, -0.00044533357857117, +0.00026500264660073, +0.00091337095312089
        };

        double[] B =
        {
            /*[ 0]*/ 0.10613121177849884,
            /*[ 1]*/ -0.7938791308754953,
            /*[ 2]*/ 2.9741147335779883,
            /*[ 3]*/ -7.369174808371755,
            /*[ 4]*/ 13.437301081373876,
            /*[ 5]*/ -18.967481201975538,
            /*[ 6]*/ 21.2267482527737,
            /*[ 7]*/ -18.967481201975538,
            /*[ 8]*/ 13.437301081373876,
            /*[ 9]*/ -7.369174808371757,
            /*[10]*/ 2.9741147335779887,
            /*[11]*/ -0.7938791308754953,
            /*[12]*/ 0.10613121177849884,
        };

        double[] A =
        {
            /*[ 0]*/ 1,
            /*[ 1]*/ -5.020620415717195,
            /*[ 2]*/ 11.694881589513276,
            /*[ 3]*/ -17.746887514727327,
            /*[ 4]*/ 20.75141302204056,
            /*[ 5]*/ -19.854306107028577,
            /*[ 6]*/ 15.390215910405962,
            /*[ 7]*/ -9.590110281996715,
            /*[ 8]*/ 4.845454719738335,
            /*[ 9]*/ -1.9390067913582834,
            /*[10]*/ 0.5682507052659562,
            /*[11]*/ -0.11013917161753861,
            /*[12]*/ 0.011626359270396867,
        };

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();
        error2.ToDebug();

        Assert.That.Value(error2).LessThan(9e-19);
    }

    [TestMethod]
    public void TypeIICorrected_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.21872794262544180, -0.43529965758160744, +0.47794280126118927, -0.00606319202674707, -0.15379106496988865,
            +0.12379156377696908, +0.22015498158263513, +0.05805388037524137, -0.01317517407603737, +0.10316624471198102,
            +0.18190835315302972, +0.10910770769374978, +0.02617691582102866, +0.05246414205404359, +0.10976184911222090,
            +0.08545437568445426, +0.01079709054953283, -0.01525428797296523, +0.01293644204213635, +0.01961964289131553,
            -0.02192353135650744, -0.05846581567937981, -0.05037273850368849, -0.02618984083737895, -0.02654782842214906,
            -0.04278595602410641, -0.04084286075696855, -0.01718276977276484, +0.00163203812124640, +0.00293579370616251,
            +0.00214194083516810, +0.01308320670784109, +0.02728725468400789, +0.03043350674727439, +0.02399132626048936,
            +0.01949009721509032, +0.01975565365483066, +0.01721975272214651, +0.00795899406342230, -0.00247800780982740,
            -0.00822454642931632, -0.01073985306018985, -0.01425604612088421, -0.01840034780721565, -0.01939182607143211,
            -0.01621539261707294, -0.01169607873921502, -0.00799194988769962, -0.00422201067559884, +0.00086548573847458,
            +0.00625082927830924, +0.00994168989905544, +0.01152799302080699, +0.01205211311596832, +0.01199908793895498,
            +0.01074536746845074, +0.00792988453627464, +0.00429570898228157, +0.00084424048862899, -0.00215579075635981,
            -0.00488970283266685, -0.00718435960755887, -0.00850181073511015, -0.00859888714812988, -0.00775882264442325,
            -0.00634694213850008, -0.00446530049280120, -0.00215622256124964, +0.00029920222466512, +0.00248283304587206,
            +0.00414431543539553, +0.00525410253824850, +0.00580688109073272, +0.00572950106078706, +0.00500678208485161,
            +0.00378636996747525, +0.00229774847217646, +0.00072047589581522, -0.00082251361106530, -0.00219146915822066,
            -0.00322534438196639, -0.00381753322217292, -0.00395738234686240, -0.00369042704721540, -0.00307130662033677,
            -0.00217102615505766, -0.00110237189876796, -0.00000535082642868, +0.00099601421593629, +0.00181439077399378,
            +0.00238765648507753, +0.00267057692704002, +0.00264785028119208, +0.00234653670839890, +0.00182669875760267,
            +0.00116121899312410, +0.00042601218797853, -0.00029966100234095, -0.00093802787217146, -0.00142620135631945
        };

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();
        error2.ToDebug();

        Assert.That.Value(error2).LessThan(3.15e-18);
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
