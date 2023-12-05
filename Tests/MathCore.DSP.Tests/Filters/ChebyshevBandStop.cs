using System.Diagnostics;

using MathCore.DSP.Filters;

using static System.Math;
using static MathCore.Polynom.Array;

using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Infrastructure;
// ReSharper disable RedundantArgumentDefaultValue
// ReSharper disable InconsistentNaming
// ReSharper disable HeuristicUnreachableCode

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ChebyshevBandStop : ChebyshevFiltersTests
{
    [TestMethod]
    public void TypeI_Even_Creation()
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
        //const double f1 = 1 / Consts.pi2;   // 0.159
        //var w1 = f1 * Consts.pi2;

        f0.AssertEquals(0.10790165633348836);
        w0.AssertEquals(0.67796610169491522);

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fpl = DigitalFilter.ToDigitalFrequency(fpl, dt);
        var Fsl = DigitalFilter.ToDigitalFrequency(fsl, dt);
        var Fsh = DigitalFilter.ToDigitalFrequency(fsh, dt);
        var Fph = DigitalFilter.ToDigitalFrequency(fph, dt);

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
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
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

        poles.AssertEquals(AccuracyComplex.Eps(1e-15),
            (-0.038245020148109342, 0.611006849626110071),
            (-0.038245020148109342, -0.611006849626110071),
            (-0.104487338181130124, 0.447288057698909503),
            (-0.104487338181130124, -0.447288057698909503),
            (-0.142732358329239445, 0.163718791927200735),
            (-0.142732358329239445, -0.163718791927200735));

        var sqrtWc = Wc.Sqrt();
        var pzf_zeros = Enumerable.Range(0, 2 * N).Select(i => i % 2 == 0 ? Complex.ImValue(sqrtWc) : Complex.ImValue(-sqrtWc));
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh).ToArray();

        //pzf_poles.ToDebugEnum();

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

        pzf_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
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

        //z_zeros.ToDebugEnum();
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

        //z_poles.ToDebugEnum()
        z_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
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

        var g_norm = (N.IsOdd() ? 1 : Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        g_norm.AssertEquals(0.0342439351432325, 1e-15);

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var h_f00 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0, dt);
        var h_fpl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fpl, dt);
        var h_fsl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsl, dt);
        var h_fcc = DoubleArrayDSPExtensions.FrequencyResponse(A, B, (fsl * fsh).Sqrt(), dt);
        var h_fsh = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsh, dt);
        var h_fph = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fph, dt);
        var h_fd5 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fd / 2, dt);

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.I);

        Assert.That.Collection(filter.A).IsEqualTo(A, 7.2e-15);
        Assert.That.Collection(filter.B).IsEqualTo(B, 5.6e-14);
    }

    [TestMethod]
    public void TypeI_Odd_Creation()
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
        const double fph = 14 / Consts.pi2; // верхняя частота границы полосы пропускания

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
        //const double f1 = 1 / Consts.pi2;   // 0.159
        //var w1 = f1 * Consts.pi2;

        f0.AssertEquals(0.12044157855602888);
        w0.AssertEquals(0.7567567567567566);

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fpl = DigitalFilter.ToDigitalFrequency(fpl, dt);
        var Fsl = DigitalFilter.ToDigitalFrequency(fsl, dt);
        var Fsh = DigitalFilter.ToDigitalFrequency(fsh, dt);
        var Fph = DigitalFilter.ToDigitalFrequency(fph, dt);

        Fpl.AssertEquals(0.31937518051807723);
        Fsl.AssertEquals(0.64524608331077715);
        Fsh.AssertEquals(2.1776750959738593);
        Fph.AssertEquals(2.681087185191323);

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
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F0 = W0 / Consts.pi2;
        const double F1 = 1 / Consts.pi2;

        W0.AssertEquals(0.7104461884009043);
        F0.AssertEquals(0.11307102268479989);

        var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.50884713990958752);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.51284645671973);
        kW.AssertEquals(1.4075661412876785);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(7);

        var beta = arcsh(1 / eps_p) / N;
        beta.AssertEquals(0.2039964798376465);

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
        //poles.ToDebugEnum();

        poles.AssertEquals(AccuracyComplex.Eps(1e-15),
            /*[ 0]*/ -0.145935804681651682,
            /*[ 1]*/ (-0.032473771555427411, 0.707095694170392974),
            /*[ 2]*/ (-0.032473771555427411, -0.707095694170392974),
            /*[ 3]*/ (-0.090989485945057832, 0.567046719980355274),
            /*[ 4]*/ (-0.090989485945057832, -0.567046719980355274),
            /*[ 5]*/ (-0.131483616730456276, 0.314687188526568629),
            /*[ 6]*/ (-0.131483616730456276, -0.314687188526568629));

        var sqrtWc = Wc.Sqrt();
        var pzf_zeros = Enumerable.Range(0, 2 * N).Select(i => i % 2 == 0 ? Complex.ImValue(sqrtWc) : Complex.ImValue(-sqrtWc));
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh).ToArray();

        //pzf_zeros.ToDebugEnum();
        pzf_zeros.AssertEquals(
            /*[ 0]*/ (0, 7.447990244669614235),
            /*[ 1]*/ (0, -7.447990244669614235),
            /*[ 2]*/ (0, 7.447990244669614235),
            /*[ 3]*/ (0, -7.447990244669614235),
            /*[ 4]*/ (0, 7.447990244669614235),
            /*[ 5]*/ (0, -7.447990244669614235),
            /*[ 6]*/ (0, 7.447990244669614235),
            /*[ 7]*/ (0, -7.447990244669614235),
            /*[ 8]*/ (0, 7.447990244669614235),
            /*[ 9]*/ (0, -7.447990244669614235),
            /*[10]*/ (0, 7.447990244669614235),
            /*[11]*/ (0, -7.447990244669614235),
            /*[12]*/ (0, 7.447990244669614235),
            /*[13]*/ (0, -7.447990244669614235));

        //pzf_poles.ToDebugEnum();
        pzf_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/ (-0.851771390245502857, 0.000000000000000000),
            /*[ 1]*/ (-65.126111677342635176, -0.000000000000000000),
            /*[ 2]*/ (-0.101685623751992776, 3.284526601999675499),
            /*[ 3]*/ (-0.522367482115460247, -16.872885543903571204),
            /*[ 4]*/ (-0.101685623751992776, -3.284526601999675499),
            /*[ 5]*/ (-0.522367482115460247, 16.872885543903571204),
            /*[ 6]*/ (-0.337732715263581618, 2.822481259261426345),
            /*[ 7]*/ (-2.318542378831646822, -19.376394756274880393),
            /*[ 8]*/ (-0.337732715263581618, -2.822481259261426345),
            /*[ 9]*/ (-2.318542378831646822, 19.376394756274880393),
            /*[10]*/ (-0.646801413466409514, 1.756830169539483677),
            /*[11]*/ (-10.237298907700751016, -27.806364057329506778),
            /*[12]*/ (-0.646801413466409514, -1.756830169539483677),
            /*[13]*/ (-10.237298907700751016, 27.806364057329506778));

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (0.756417559622531210, 0.654089042481751148),
            /*[ 1]*/ (0.756417559622531210, -0.654089042481751148),
            /*[ 2]*/ (0.756417559622531210, 0.654089042481751148),
            /*[ 3]*/ (0.756417559622531210, -0.654089042481751148),
            /*[ 4]*/ (0.756417559622531210, 0.654089042481751148),
            /*[ 5]*/ (0.756417559622531210, -0.654089042481751148),
            /*[ 6]*/ (0.756417559622531210, 0.654089042481751148),
            /*[ 7]*/ (0.756417559622531210, -0.654089042481751148),
            /*[ 8]*/ (0.756417559622531210, 0.654089042481751148),
            /*[ 9]*/ (0.756417559622531210, -0.654089042481751148),
            /*[10]*/ (0.756417559622531210, 0.654089042481751148),
            /*[11]*/ (0.756417559622531210, -0.654089042481751148),
            /*[12]*/ (0.756417559622531210, 0.654089042481751148),
            /*[13]*/ (0.756417559622531210, -0.654089042481751148));

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/ 0.918302251228021560,
            /*[ 1]*/ -0.530108926487634968,
            /*[ 2]*/ (0.938138236661072322, 0.316683223278757087),
            /*[ 3]*/ (0.162967963632972723, -0.956157975375219626),
            /*[ 4]*/ (0.938138236661072322, -0.316683223278757087),
            /*[ 5]*/ (0.162967963632972723, 0.956157975375219626),
            /*[ 6]*/ (0.929622969434593438, 0.267794092139977180),
            /*[ 7]*/ (0.021955725315704035, -0.887236147461566760),
            /*[ 8]*/ (0.929622969434593438, -0.267794092139977180),
            /*[ 9]*/ (0.021955725315704035, 0.887236147461566760),
            /*[10]*/ (0.923420034559933933, 0.163663236631297865),
            /*[11]*/ (-0.283258983144856313, -0.659118452026115409),
            /*[12]*/ (0.923420034559933933, -0.163663236631297865),
            /*[13]*/ (-0.283258983144856313, 0.659118452026115409));

        var g_norm = (N.IsOdd() ? 1 : Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        g_norm.AssertEquals(0.029318326737782795, 1e-15);

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var h_f00 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0, dt);
        var h_fpl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fpl, dt);
        var h_fsl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsl, dt);
        var h_fcc = DoubleArrayDSPExtensions.FrequencyResponse(A, B, (fsl * fsh).Sqrt(), dt);
        var h_fsh = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsh, dt);
        var h_fph = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fph, dt);
        var h_fd5 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fd / 2, dt);

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.I);

        Assert.That.Collection(filter.A).IsEqualTo(A, 2.2e-14);
        Assert.That.Collection(filter.B).IsEqualTo(B, 6.4e-14);
    }

    [TestMethod]
    public void TypeII_Even_Creation()
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
        //const double f1 = 1 / Consts.pi2;   // 0.159
        //var w1 = f1 * Consts.pi2;

        f0.AssertEquals(0.10790165633348836);
        w0.AssertEquals(0.67796610169491522);

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fpl = DigitalFilter.ToDigitalFrequency(fpl, dt);
        var Fsl = DigitalFilter.ToDigitalFrequency(fsl, dt);
        var Fsh = DigitalFilter.ToDigitalFrequency(fsh, dt);
        var Fph = DigitalFilter.ToDigitalFrequency(fph, dt);

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
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
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
        if (N.IsOdd()) poles[0] = -W0 / sh;
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

        poles.AssertEquals(AccuracyComplex.Eps(1e-15),
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

        var pzf_zeros_enum = AnalogBasedFilter.TransformToBandStop(zeros, Fsl, Fsh);
        if (N.IsOdd())
        {
            var sqrt_wc = (Wsl * Wsh).Sqrt();
            pzf_zeros_enum = pzf_zeros_enum.AppendLast(
                Complex.ImValue(+sqrt_wc),
                Complex.ImValue(-sqrt_wc)
            );
        }
        var pzf_zeros = pzf_zeros_enum.ToArray();
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

        pzf_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
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
        z_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
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

        var g_norm = 1 / zeros_to_poles;

        g_norm.AssertEquals(0.10613121177849884, 1e-15);

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        //A.ToDebugEnum();

        var h_f00 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0, dt);
        var h_fpl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fpl, dt);
        var h_fsl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsl, dt);
        var h_fcc = DoubleArrayDSPExtensions.FrequencyResponse(A, B, (fsl * fsh).Sqrt(), dt);
        var h_fsh = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsh, dt);
        var h_fph = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fph, dt);
        var h_fd5 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fd / 2, dt);

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.II);

        //B.ToDebugEnum();
        //A.ToDebugEnum();

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);
    }

    [TestMethod]
    public void TypeII_Odd_Creation()
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
        const double fph = 14 / Consts.pi2; // верхняя частота границы полосы пропускания

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
        //const double f1 = 1 / Consts.pi2;   // 0.159
        //var w1 = f1 * Consts.pi2;

        f0.AssertEquals(0.12044157855602888);
        w0.AssertEquals(0.7567567567567566);

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fpl = DigitalFilter.ToDigitalFrequency(fpl, dt);
        var Fsl = DigitalFilter.ToDigitalFrequency(fsl, dt);
        var Fsh = DigitalFilter.ToDigitalFrequency(fsh, dt);
        var Fph = DigitalFilter.ToDigitalFrequency(fph, dt);

        Fpl.AssertEquals(0.31937518051807723);
        Fsl.AssertEquals(0.64524608331077715);
        Fsh.AssertEquals(2.1776750959738593);
        Fph.AssertEquals(2.681087185191323);

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
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F0 = W0 / Consts.pi2;
        const double F1 = 1 / Consts.pi2;

        W0.AssertEquals(0.7104461884009043);
        F0.AssertEquals(0.11307102268479989);

        var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.50884713990958752);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F1 / F0;

        kEps.AssertEquals(196.51284645671973);
        kW.AssertEquals(1.4075661412876785);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(7);

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(0.7568989093729265);

        var poles = new Complex[N];
        var r = N % 2;
        var sh = Sinh(beta);
        var ch = Cosh(beta);

        if (N.IsOdd()) poles[0] = -W0 / sh;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (cos, sin) = Complex.Exp(dth * (i - r + 1));
            var z = new Complex(-sh * sin, ch * cos);
            var norm = W0 / z.Power;
            (poles[i], poles[i + 1]) = Complex.Conjugate(z.Re * norm, z.Im * norm);
        }

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ -0.854653317262274226,
            /*[ 1]*/ (-0.080057986392564895, 0.548703565842303176),
            /*[ 2]*/ (-0.080057986392564895, -0.548703565842303176),
            /*[ 3]*/ (-0.282749324575917005, 0.554647252362557519),
            /*[ 4]*/ (-0.282749324575917005, -0.554647252362557519),
            /*[ 5]*/ (-0.605151136846373894, 0.455888810284692192),
            /*[ 6]*/ (-0.605151136846373894, -0.455888810284692192)
        );

        var L = N / 2;
        var zeros = new Complex[2 * L];
        for (var (n, dth) = (1, PI / N); n <= L; n++)
            (zeros[2 * n - 2], zeros[2 * n - 1]) = Complex.Conjugate(0, W0 / Cos(dth * (n - 0.5)));

        //zeros.AsEnumerable().ToIm().ToDebugEnum();
        zeros.AsEnumerable().ToIm().AssertEquals(
            /*[ 0]*/ +0.7287166358905175,
            /*[ 1]*/ -0.7287166358905175,
            /*[ 2]*/ +0.9086947818450832,
            /*[ 3]*/ -0.9086947818450832,
            /*[ 4]*/ +1.6374114177356005,
            /*[ 5]*/ -1.6374114177356005
        );

        var pzf_zeros_enum = AnalogBasedFilter.TransformToBandStop(zeros, Fsl, Fsh);
        if (N.IsOdd())
        {
            var sqrt_wc = (Wsl * Wsh).Sqrt();
            pzf_zeros_enum = pzf_zeros_enum.AppendLast(
                Complex.ImValue(+sqrt_wc),
                Complex.ImValue(-sqrt_wc)
            );
        }
        var pzf_zeros = pzf_zeros_enum.ToArray();
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh).ToArray();

        //pzf_zeros.AsEnumerable().ToRe().Sum(v => v * v).ToDebug();
        pzf_zeros.AsEnumerable().ToRe().Sum(v => v * v).AssertLessThan(3.8e-30);

        //pzf_zeros.AsEnumerable().ToIm().ToDebugEnum();
        pzf_zeros.AsEnumerable().ToIm().AssertEquals(
            /*[ 0]*/ +03.349321876161251,
            /*[ 1]*/ -16.562325370851593,
            /*[ 2]*/ +16.562325370851593,
            /*[ 3]*/ -03.349321876161251,
            /*[ 4]*/ +03.842095143156990,
            /*[ 5]*/ -14.438101248870378,
            /*[ 6]*/ +14.438101248870378,
            /*[ 7]*/ -03.842095143156990,
            /*[ 8]*/ +05.067149123771191,
            /*[ 9]*/ -10.947488879784274,
            /*[10]*/ +10.947488879784274,
            /*[11]*/ -05.067149123771191,
            /*[12]*/ +07.447990244669614,
            /*[13]*/ -07.447990244669614
        );

        //pzf_poles.ToDebugEnum();
        pzf_poles.AssertEquals(
            /*[ 0]*/ (-5.633006543228341378, +04.872555383845363686),
            /*[ 1]*/ (-5.633006543228341378, -04.872555383845363686),
            /*[ 2]*/ (-0.303905998847857517, +02.749558700682777967),
            /*[ 3]*/ (-2.203017755114543785, -19.931579696675406410),
            /*[ 4]*/ (-0.303905998847857517, -02.749558700682777967),
            /*[ 5]*/ (-2.203017755114543785, +19.931579696675406410),
            /*[ 6]*/ (-1.051002269785872301, +02.942088876784905160),
            /*[ 7]*/ (-5.973241997226903521, -16.720999890862426440),
            /*[ 8]*/ (-1.051002269785872301, -02.942088876784905160),
            /*[ 9]*/ (-5.973241997226903521, +16.720999890862426440),
            /*[10]*/ (-2.374573487550539763, +03.361801543916515644),
            /*[11]*/ (-7.775753568728177356, -11.008520262487085262),
            /*[12]*/ (-2.374573487550539763, -03.361801543916515644),
            /*[13]*/ (-7.775753568728177356, +11.008520262487085262)
        );

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (0.945440334798850590, 0.325795293605412928),
            /*[ 1]*/ (0.186396853715641519, -0.982474535509653801),
            /*[ 2]*/ (0.186396853715641519, 0.982474535509653801),
            /*[ 3]*/ (0.945440334798850590, -0.325795293605412928),
            /*[ 4]*/ (0.928818426272944642, 0.370535195380749471),
            /*[ 5]*/ (0.314797390561100621, -0.949158892333587079),
            /*[ 6]*/ (0.314797390561100621, 0.949158892333587079),
            /*[ 7]*/ (0.928818426272944642, -0.370535195380749471),
            /*[ 8]*/ (0.879363646398171994, 0.476150792704696402),
            /*[ 9]*/ (0.538912816438409492, -0.842361547245849218),
            /*[10]*/ (0.538912816438409492, 0.842361547245849218),
            /*[11]*/ (0.879363646398171994, -0.476150792704696402),
            /*[12]*/ (0.756417559622531210, 0.654089042481751148),
            /*[13]*/ (0.756417559622531210, -0.654089042481751148)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.506067971025137342, 0.286287119237451870),
            /*[ 1]*/ (0.506067971025137342, -0.286287119237451870),
            /*[ 2]*/ (0.934586620991862072, 0.261982077550719172),
            /*[ 3]*/ (-0.002382674320049361, -0.895557957611133770),
            /*[ 4]*/ (0.934586620991862072, -0.261982077550719172),
            /*[ 5]*/ (-0.002382674320049361, 0.895557957611133770),
            /*[ 6]*/ (0.863742813761753880, 0.260476766440063345),
            /*[ 7]*/ (0.088795264844538230, -0.700942358546569944),
            /*[ 8]*/ (0.863742813761753880, -0.260476766440063345),
            /*[ 9]*/ (0.088795264844538230, 0.700942358546569944),
            /*[10]*/ (0.748275629893491945, 0.262680122820566186),
            /*[11]*/ (0.244600682130233199, -0.493279571840732067),
            /*[12]*/ (0.748275629893491945, -0.262680122820566186),
            /*[13]*/ (0.244600682130233199, 0.493279571840732067)
        );

        var zeros_to_poles = z_zeros.Multiply(z => 1 - z).Re / z_poles.Multiply(z => 1 - z).Re;

        var G_norm = 1 / zeros_to_poles;

        G_norm.AssertEquals(0.14070756712838955);

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        //A.ToDebugEnum();

        var h_f00 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0, dt);
        var h_fpl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fpl, dt);
        var h_fsl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsl, dt);
        var h_fcc = DoubleArrayDSPExtensions.FrequencyResponse(A, B, (fsl * fsh).Sqrt(), dt);
        var h_fsh = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsh, dt);
        var h_fph = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fph, dt);
        var h_fd5 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fd / 2, dt);

        //var h_f00_db = h_f00.Power.In_dB_byPower().ToDebug();
        //var h_fpl_db = h_fpl.Power.In_dB_byPower().ToDebug();
        //var h_fsl_db = h_fsl.Power.In_dB_byPower().ToDebug();
        //var h_fcc_db = h_fcc.Power.In_dB_byPower().ToDebug();
        //var h_fsh_db = h_fsh.Power.In_dB_byPower().ToDebug();
        //var h_fph_db = h_fph.Power.In_dB_byPower().ToDebug();
        //var h_fd5_db = h_fd5.Power.In_dB_byPower().ToDebug();

        h_f00.Abs.AssertGreaterOrEqualsThan(Gp, 3.4e-11);
        h_fpl.Abs.AssertGreaterOrEqualsThan(Gp);

        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fcc.Abs.AssertLessOrEqualsThan(Gs);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);

        h_fph.Abs.AssertLessOrEqualsThan(Gp);
        h_fd5.Abs.AssertGreaterOrEqualsThan(Gp, 1.9e-15);

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.II);

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);
    }

    [TestMethod]
    public void TypeII_Corrected_Even_Creation()
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
        //const double f1 = 1 / Consts.pi2;   // 0.159
        //var w1 = f1 * Consts.pi2;

        f0.AssertEquals(0.10790165633348836);
        w0.AssertEquals(0.67796610169491522);

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fpl = DigitalFilter.ToDigitalFrequency(fpl, dt);
        var Fsl = DigitalFilter.ToDigitalFrequency(fsl, dt);
        var Fsh = DigitalFilter.ToDigitalFrequency(fsh, dt);
        var Fph = DigitalFilter.ToDigitalFrequency(fph, dt);

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
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
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
        poles.AssertEquals(AccuracyComplex.Eps(1e-15),
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

        var h_f00 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0, dt);
        var h_fpl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fpl, dt);
        var h_fsl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsl, dt);
        var h_fcc = DoubleArrayDSPExtensions.FrequencyResponse(A, B, (fsl * fsh).Sqrt(), dt);
        var h_fsh = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsh, dt);
        var h_fph = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fph, dt);
        var h_fd5 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fd / 2, dt);

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.IICorrected);

        filter.A.AssertEquals(A);
        filter.B.AssertEquals(B);
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_Creation()
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
        const double fph = 14 / Consts.pi2; // верхняя частота границы полосы пропускания

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
        //const double f1 = 1 / Consts.pi2;   // 0.159
        //var w1 = f1 * Consts.pi2;

        f0.AssertEquals(0.12044157855602888);
        w0.AssertEquals(0.7567567567567566);

        // Преобразуем частоты аналогового фильтра в частоты цифрового фильтра с учётом заданной частоты дискретизации
        var Fpl = DigitalFilter.ToDigitalFrequency(fpl, dt);
        var Fsl = DigitalFilter.ToDigitalFrequency(fsl, dt);
        var Fsh = DigitalFilter.ToDigitalFrequency(fsh, dt);
        var Fph = DigitalFilter.ToDigitalFrequency(fph, dt);

        Fpl.AssertEquals(0.31937518051807723);
        Fsl.AssertEquals(0.64524608331077715);
        Fsh.AssertEquals(2.1776750959738593);
        Fph.AssertEquals(2.681087185191323);

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
        //const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F0 = W0 / Consts.pi2;
        const double F1 = 1 / Consts.pi2;

        W0.AssertEquals(0.7104461884009043);
        F0.AssertEquals(0.11307102268479989);

        var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.50884713990958752);
        eps_s.AssertEquals(99.994999874993752);

        var kEps = eps_s / eps_p;
        var kW = F0 / F1;

        kEps.AssertEquals(196.51284645671973);
        kW.AssertEquals(0.7104461884009043);

        var N = (int)Ceiling(arcch(kEps) / arcch(1 / kW)); // Порядок фильтра
        N.AssertEquals(7);
        //N.ToDebug();

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(0.7568989093729265);
        //beta.ToDebug();

        var r = N % 2;
        var sh = Sinh(beta);
        var ch = Cosh(beta);

        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -1 / sh;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (sin, cos) = Complex.SinCos(dth * (i - r + 1));
            var z = new Complex(-sin * sh, cos * ch);
            var norm = 1 / z.Power;
            (poles[i], poles[i + 1]) = Complex.Conjugate(z.Re * norm, z.Im * norm);
        }

        var zeros = new Complex[N - r];
        for (var (n, dth, L) = (1, PI / N, N / 2); n <= L; n++)
        {
            var th = dth * (n - 0.5);
            (zeros[2 * n - 2], zeros[2 * n - 1]) = Complex.Conjugate(0, W0 / Cos(th) / kW);
        }

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ -01.202981071917573530,
            /*[ 1]*/ (-0.112686910985844033, +0.772336560883440271),
            /*[ 2]*/ (-0.112686910985844033, -0.772336560883440271),
            /*[ 3]*/ (-0.397988375745020939, +0.780702692783778351),
            /*[ 4]*/ (-0.397988375745020939, -0.780702692783778351),
            /*[ 5]*/ (-0.851790250586702435, +0.641693653548654708),
            /*[ 6]*/ (-0.851790250586702435, -0.641693653548654708)
        );

        //poles.Zip(poles1, (p, p1) => p - p1).ToDebugEnum();

        //zeros.AsEnumerable().ToRe().Sum(v => v * v).ToDebug();
        zeros.AsEnumerable().ToRe().Sum(v => v * v).AssertEquals(0);

        //zeros.AsEnumerable().ToIm().ToDebugEnum();
        zeros.AsEnumerable().ToIm().AssertEquals(
            /*[ 0]*/ +1.0257168632725540,
            /*[ 1]*/ -1.0257168632725540,
            /*[ 2]*/ +1.2790480076899327,
            /*[ 3]*/ -1.2790480076899327,
            /*[ 4]*/ +2.3047648709624860,
            /*[ 5]*/ -2.3047648709624860
        );

        var pzf_zeros_enum = AnalogBasedFilter.TransformToBandStop(zeros, Fsl, Fsh);
        if (N.IsOdd())
        {
            var sqrt_wc = (Wsl * Wsh).Sqrt();
            pzf_zeros_enum = pzf_zeros_enum.AppendLast(
                Complex.ImValue(+sqrt_wc),
                Complex.ImValue(-sqrt_wc)
            );
        }
        var pzf_zeros = pzf_zeros_enum.ToArray();
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh).ToArray();

        //pzf_zeros.AsEnumerable().ToRe().Sum(v => v * v).ToDebug();
        pzf_zeros.AsEnumerable().ToRe().Sum(v => v * v).AssertEquals(0, 3.2e-30);

        //pzf_zeros.AsEnumerable().ToIm().ToDebugEnum();
        pzf_zeros.AsEnumerable().ToIm().AssertEquals(
            /*[ 0]*/ +04.109963805916897,
            /*[ 1]*/ -13.497091776047480,
            /*[ 2]*/ +13.497091776047480,
            /*[ 3]*/ -04.109963805916897,
            /*[ 4]*/ +04.581103280147630,
            /*[ 5]*/ -12.108995430224418,
            /*[ 6]*/ +12.108995430224418,
            /*[ 7]*/ -04.581103280147630,
            /*[ 8]*/ +05.646526579500983,
            /*[ 9]*/ -09.824191545662782,
            /*[10]*/ +09.824191545662782,
            /*[11]*/ -05.646526579500983,
            /*[12]*/ +07.447990244669614,
            /*[13]*/ -07.447990244669614
        );

        //pzf_poles.ToDebugEnum();
        pzf_poles.AssertEquals(
            /*[ 0]*/ (-4.001948027873929448, +06.281478382267164484),
            /*[ 1]*/ (-4.001948027873929448, -06.281478382267164484),
            /*[ 2]*/ (-0.324629645893097840, +03.501333380056669498),
            /*[ 3]*/ (-1.456404779721176546, -15.708234705683944199),
            /*[ 4]*/ (-0.324629645893097840, -03.501333380056669498),
            /*[ 5]*/ (-1.456404779721176546, +15.708234705683944199),
            /*[ 6]*/ (-1.086714410588036417, +03.776480011722634700),
            /*[ 7]*/ (-3.903633155308095404, -13.565654821989250323),
            /*[ 8]*/ (-1.086714410588036417, -3.7764800117226347000),
            /*[ 9]*/ (-3.903633155308095404, +13.565654821989250323),
            /*[10]*/ (-2.243591127622733161, +04.474353915760326927),
            /*[11]*/ (-4.967670040533056053, -09.906936083142635852),
            /*[12]*/ (-2.243591127622733161, -04.474353915760326927),
            /*[13]*/ (-4.967670040533056053, +09.906936083142635852)
        );

        //pzf_poles.Zip(pzf_poles1, (p, p1) => p - p1).ToDebugEnum();

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (0.918963134390169301, +0.394343451361621911),
            /*[ 1]*/ (0.374165548594304176, -0.927361926242997736),
            /*[ 2]*/ (0.374165548594304176, +0.927361926242997736),
            /*[ 3]*/ (0.918963134390169301, -0.394343451361621911),
            /*[ 4]*/ (0.900298432841460605, +0.435273169197470833),
            /*[ 5]*/ (0.463518766194801424, -0.886087102595026432),
            /*[ 6]*/ (0.463518766194801424, +0.886087102595026432),
            /*[ 7]*/ (0.900298432841460605, -0.435273169197470833),
            /*[ 8]*/ (0.852352405439963490, +0.522967854595966841),
            /*[ 9]*/ (0.611231120276368167, -0.791452157496393416),
            /*[10]*/ (0.611231120276368167, +0.791452157496393416),
            /*[11]*/ (0.852352405439963490, -0.522967854595966841),
            /*[12]*/ (0.756417559622531210, +0.654089042481751148),
            /*[13]*/ (0.756417559622531210, -0.654089042481751148)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.559706234698156502, +0.408186076587035207),
            /*[ 1]*/ (0.559706234698156502, -0.408186076587035207),
            /*[ 2]*/ (0.911332705091532747, +0.329266171995574786),
            /*[ 3]*/ (0.213724653872510523, -0.888567862455853397),
            /*[ 4]*/ (0.911332705091532747, -0.329266171995574786),
            /*[ 5]*/ (0.213724653872510523, +0.888567862455853397),
            /*[ 6]*/ (0.837977249191020901, +0.329168603909464119),
            /*[ 7]*/ (0.265729179818831762, -0.718319472189897001),
            /*[ 8]*/ (0.837977249191020901, -0.329168603909464119),
            /*[ 9]*/ (0.265729179818831762, +0.718319472189897001),
            /*[10]*/ (0.728338195589888215, +0.347659545116905500),
            /*[11]*/ (0.384147756131204554, -0.549216780234427437),
            /*[12]*/ (0.728338195589888215, -0.347659545116905500),
            /*[13]*/ (0.384147756131204554, +0.549216780234427437)
        );

        var re_to_im = z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z);

        var G_norm = 1 / re_to_im.Re;
        //G_norm.ToDebug();

        G_norm.AssertEquals(0.22881500023034296);

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b.Re * G_norm);
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var h_f00 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0, dt);
        var h_fpl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fpl, dt);
        var h_fsl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsl, dt);
        var h_fcc = DoubleArrayDSPExtensions.FrequencyResponse(A, B, (fsl * fsh).Sqrt(), dt);
        var h_fsh = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsh, dt);
        var h_fph = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fph, dt);
        var h_fd5 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fd / 2, dt);

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.IICorrected);

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);
    }

    [TestMethod]
    public void TypeI_Even_ImpulseResponse()
    {
        double[] expected_h =
        [
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
        ];

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.I);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();
        error2.ToDebug();

        Assert.That.Value(error2).LessThan(2e-19);
    }

    [TestMethod]
    public void TypeI_Odd_ImpulseResponse()
    {
        double[] expected_h =
        [
            +0.02931832673469236, -0.14119590714658910, +0.33769354399014430, -0.42152072005675090, +0.22647622495935696,
            +0.12504817090184703, -0.18003417501798627, -0.03309124808808608, +0.19127649378592768, +0.02082762641067514,
            -0.08714060810652260, +0.03013115279466035, +0.13317444696585004, +0.07780377522829712, -0.01074200127329353,
            +0.02704355226816512, +0.13345555699567463, +0.12866413352480308, +0.03740788024813411, +0.03297921729155228,
            +0.11454498379981340, +0.13050926731209353, +0.05825338506239403, +0.01891510109260321, +0.05582526672615626,
            +0.08066413269873648, +0.03672565933969205, -0.01773168981468920, -0.01845153195815877, +0.00639444539394439,
            -0.00405764981869209, -0.04486795185113604, -0.06255968443795577, -0.04169929004361913, -0.02190913634555257,
            -0.03191960139185412, -0.04833188636435190, -0.03615528697706023, -0.00354102030518446, +0.01218843472512218,
            +0.00263011481980758, -0.00074678459422150, +0.02111123125480063, +0.04357493760788578, +0.03825037041352632,
            +0.01866300615914991, +0.01712012703964886, +0.03164567923497831, +0.03070708811568959, +0.00664621072832370,
            -0.01156558845372643, -0.00613145609363948, +0.00192836465478171, -0.01009300472217869, -0.02854053229199179,
            -0.02698658575852132, -0.00981916644708191, -0.00348937941721065, -0.01255051723615001, -0.01413075957949138,
            +0.00248207538874954, +0.01829508297700273, +0.01631488971962101, +0.00810368621705589, +0.01259399384706852,
            +0.02454899660677456, +0.02433659727513267, +0.01057803824684508, +0.00150604001867849, +0.00482859112042983,
            +0.00664050611476313, -0.00404123919272170, -0.01698389842973665, -0.01783255493604002, -0.01090864020393154,
            -0.01069042484653378, -0.01766474802503874, -0.01842346162917170, -0.00836594794359846, +0.00097796555476736,
            +0.00093694138056376, -0.00112881399342169, +0.00469281893325909, +0.01411759436576898, +0.01592285816788808,
            +0.01013335495142117, +0.00716538185476216, +0.01043229320153959, +0.01145452180996506, +0.00461295780096530,
            -0.00347211128922250, -0.00482348394339025, -0.00254277827741526, -0.00467284619984420, -0.01053575706245073,
            -0.01209264478285616, -0.00726618237737667, -0.00300380562597034, -0.00372453913143892, -0.00464703830562346
        ];

        const double fd = 10;         // Частота дискретизации
        const double dt = 1 / fd;       // Период дискретизации

        const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
        const double Rs = 40;   // Уровень подавления более 40 дБ

        var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
        var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

        const double fpl = 2 / Consts.pi2;  // нижняя частота границы полосы пропускания
        const double fsl = 4 / Consts.pi2;  // нижняя частота границы полосы заграждения
        const double fsh = 12 / Consts.pi2; // верхняя частота границы полосы заграждения
        const double fph = 14 / Consts.pi2; // верхняя частота границы полосы пропускания

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.I);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();
        error2.ToDebug();

        Assert.That.Value(error2).LessThan(1.21e-17);
    }

    [TestMethod]
    public void TypeII_Even_ImpulseResponse()
    {
        double[] expected_h =
        [
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
        ];

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.II);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();
        error2.ToDebug();

        Assert.That.Value(error2).LessThan(9.01e-19);
    }

    [TestMethod]
    public void TypeII_Odd_ImpulseResponse()
    {
        double[] expected_h =
        [
            +0.14070756712330665, -0.32825969995649673, +0.47106622467692844, -0.23043531279864370, -0.03854744867834650,
            +0.23182995956020150, +0.02991639158078585, -0.07252318793021478, +0.10797864809945124, +0.15869538680745160,
            +0.02493429637195324, +0.01415841365118995, +0.12724577678351554, +0.12817037135712680, +0.03813600964044966,
            +0.04118796331741659, +0.10718030959279506, +0.09021392671061740, +0.02292884611175169, +0.02343697352817929,
            +0.06152791760342690, +0.04162193102694615, -0.00945462066028003, -0.01337432980780893, +0.00835060028738142,
            -0.00663828390770037, -0.03978064602489056, -0.04017229813325618, -0.02251912180897142, -0.02750929506095161,
            -0.04315612236878769, -0.03699675600909061, -0.01914720381450367, -0.01614237691259401, -0.02037867639741215,
            -0.01128817642517591, +0.00421484916113197, +0.00893359194710357, +0.00779777423602078, +0.01395777137136056,
            +0.02288281467881818, +0.02365819001137046, +0.01972732458481696, +0.01974699059848409, +0.02099681273862940,
            +0.01678327061681784, +0.00969276975246535, +0.00556445645112601, +0.00285586892250456, -0.00251160193990123,
            -0.00868787140500440, -0.01186069877782131, -0.01306212965070205, -0.01498482673548182, -0.01655978566555560,
            -0.01554703495877278, -0.01287297182793961, -0.01052236316079329, -0.00805082251799172, -0.00426549185782252,
            +0.00001796272599173, +0.00337141819734970, +0.00599397826500036, +0.00861431606260068, +0.01072754333765069,
            +0.01150524773129065, +0.01121857817644415, +0.01050634849433266, +0.00923339385921275, +0.00707138130156770,
            +0.00439341967168809, +0.00178076655209992, -0.00071851420765528, -0.00321747362609982, -0.00542341384399943,
            -0.00696194417143695, -0.00784847155984053, -0.00823878146115010, -0.00805850451019054, -0.00720491019589881,
            -0.00583701559035437, -0.00420093398386269, -0.00238242059785154, -0.00042970744235579, +0.00146833665995138,
            +0.00310107163760347, +0.00439626720179740, +0.00533900170322450, +0.00585477125153213, +0.00588187297455566,
            +0.00546692907072331, +0.00470394244771775, +0.00365141939072080, +0.00237168355087302, +0.00098471595732180,
            -0.00037750530157191, -0.00162817418287264, -0.00270224227578939, -0.00352248632076450, -0.00402781595501177
        ];

        const double fd = 10;           // Частота дискретизации
        const double dt = 1 / fd;       // Период дискретизации

        const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
        const double Rs = 40;   // Уровень подавления более 40 дБ

        var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
        var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

        const double fpl = 2 / Consts.pi2;  // нижняя частота границы полосы пропускания
        const double fsl = 4 / Consts.pi2;  // нижняя частота границы полосы заграждения
        const double fsh = 12 / Consts.pi2; // верхняя частота границы полосы заграждения
        const double fph = 14 / Consts.pi2; // верхняя частота границы полосы пропускания

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.II);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average().Sqrt();
        error.ToDebug();

        Assert.That.Value(error).LessThan(5.36e-12);
    }

    [TestMethod]
    public void TypeII_Corrected_Even_ImpulseResponse()
    {
        double[] expected_h =
        [
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
        ];

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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.IICorrected);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();
        error2.ToDebug();

        Assert.That.Value(error2).LessThan(3.15e-18);
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_ImpulseResponse()
    {
        double[] expected_h =
        [
            +0.22881500022931373, -0.44664275856193460, +0.47358299185039120, +0.02307987508704328, -0.15678168951627006,
            +0.10218005448231729, +0.23153945606384074, +0.08734899468490825, -0.01440115295783170, +0.08020284319872750,
            +0.18380972040665028, +0.13544097316529410, +0.03185074186421728, +0.02701393820572644, +0.09347467595244445,
            +0.09919652249784168, +0.02220189399682204, -0.03554405564211223, -0.01713849711003434, +0.01536584487112376,
            -0.00774868836672680, -0.06010997069818232, -0.07177526859472418, -0.03685495602038790, -0.01236147966455513,
            -0.02698923454290885, -0.04473288188041932, -0.02700884161234825, +0.01060802412713162, +0.02714734080237191,
            +0.01661361399901176, +0.00966779653657926, +0.02361645056853234, +0.03961968939724947, +0.03493561758837013,
            +0.01585963553576044, +0.00477255452260211, +0.00690332363299641, +0.00681140883365204, -0.00541263184543605,
            -0.02006876855209103, -0.02381237262454140, -0.01834517891928529, -0.01514153339479393, -0.01760459374985310,
            -0.01783167256219973, -0.01010270768944888, +0.00070638555392700, +0.00687276872494584, +0.00826486904657544,
            +0.01043482685867396, +0.01517629802500484, +0.01815170389947627, +0.01580771680094576, +0.01032618355055682,
            +0.00594563800008821, +0.00320519343540165, -0.00051540729980875, -0.00604265571242141, -0.01083109815374692,
            -0.01250596706202954, -0.01175919221920397, -0.01066856612993569, -0.00957336331041529, -0.00712798763737360,
            -0.00289512944074976, +0.00164863981692409, +0.00491787818256083, +0.00690093733314759, +0.00843657374403079,
            +0.00955071587439764, +0.00938215599752991, +0.00755065630171347, +0.00479265898352213, +0.00205986801628473,
            -0.00044215488484536, -0.00296571701584633, -0.00535019788396005, -0.00694116704651249, -0.00732316431107390,
            -0.00672590698429617, -0.00559455692367545, -0.00405862765446835, -0.00204091065494762, +0.00028665487048555,
            +0.00244000298377643, +0.00401910824644050, +0.00496866090419134, +0.00538941174728274, +0.00525956689912978,
            +0.00448032749233809, +0.00312234928628582, +0.00147601793156162, -0.00015860892867822, -0.00163441644596733,
            -0.00288335645194253, -0.00377852000072013, -0.00416048380510373, -0.00397910346586083, -0.00333789128109529,
        ];

        const double fd = 10;           // Частота дискретизации
        const double dt = 1 / fd;       // Период дискретизации

        const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
        const double Rs = 40;   // Уровень подавления более 40 дБ

        var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
        var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

        const double fpl = 2 / Consts.pi2;  // нижняя частота границы полосы пропускания
        const double fsl = 4 / Consts.pi2;  // нижняя частота границы полосы заграждения
        const double fsh = 12 / Consts.pi2; // верхняя частота границы полосы заграждения
        const double fph = 14 / Consts.pi2; // верхняя частота границы полосы пропускания

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.IICorrected);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average().Sqrt();
        //error2.ToDebug();

        Assert.That.Value(error2).LessThan(5.9e-10);
    }

    [TestMethod]
    public void TypeI_Even_SignalProcessing()
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

        // Расчёт коэффициентов передачи по мощности
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

        // Расчёт коэффициентов передачи по мощности в логарифмическим масштабе
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
        // должны быть ниже уровня подавления Rs
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

    [TestMethod]
    public void TypeI_Odd_SignalProcessing()
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
        const double fph = 14 / Consts.pi2; // верхняя частота границы полосы пропускания

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

        // Расчёт коэффициентов передачи по мощности
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

        // Расчёт коэффициентов передачи по мощности в логарифмическим масштабе
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
        // должны быть ниже уровня подавления Rs
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

    [TestMethod]
    public void TypeII_Even_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.II);

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

        // Расчёт коэффициентов передачи по мощности
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

        // Расчёт коэффициентов передачи по мощности в логарифмическим масштабе
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
        // должны быть ниже уровня подавления Rs
        k_sl_db.AssertLessOrEqualsThan(-Rs, 2.16);      // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_sl_sh_db.AssertLessOrEqualsThan(-Rs);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertLessOrEqualsThan(-Rs, 1.5);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_sh_sl_db.AssertLessOrEqualsThan(-Rs, 1.5);    // Коэффициент передачи у верхнего края полосы подавления
        k_sh_db.AssertLessOrEqualsThan(-Rs, 0.36);      // Коэффициент передачи на верхней границе полосы подавления

        // Коэффициенты передачи в переходной полосе
        // между верхней полосой заграждения и верхней
        // полосой пропускания
        k_sh_ph_db.AssertLessThan(-Rp);
        k_ph_sh_db.AssertLessThan(-Rs);
        k_ph_sh_db.AssertLessThan(-Rp);

        // Коэффициенты передачи в верхней полосе пропускания
        // от верхней границы полосы пропускания до частоты,
        // близкой к половине частоты дискретизации
        k_ph_db.AssertLessThan(-Rp, 0.01);
        k_ph_fd05_db.AssertGreaterOrEqualsThan(-Rs);
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
        h_sl_pl.AssertLessOrEqualsThan(-Rp);

        h_sl.AssertLessOrEqualsThan(-Rs);
        h_sl_sh.AssertLessOrEqualsThan(-Rs);
        h_c0.AssertLessOrEqualsThan(-Rs);
        h_sh_sl.AssertLessOrEqualsThan(-Rs);
        h_sh.AssertLessOrEqualsThan(-Rs);

        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_ph_sh.AssertLessOrEqualsThan(-Rp);

        h_ph.AssertLessOrEqualsThan(-Rp);
        h_ph_fd05.AssertLessOrEqualsThan(-Rp);
        h_fd05.AssertGreaterOrEqualsThan(-Rp);
    }

    [TestMethod]
    public void TypeII_Odd_SignalProcessing()
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
        const double fph = 14 / Consts.pi2; // верхняя частота границы полосы пропускания

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.II);

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

        // Расчёт коэффициентов передачи по мощности
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

        // Расчёт коэффициентов передачи по мощности в логарифмическим масштабе
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
        // должны быть ниже уровня подавления Rs
        k_sl_db.AssertLessOrEqualsThan(-Rs, 2.16);      // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_sl_sh_db.AssertLessOrEqualsThan(-Rs);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertLessOrEqualsThan(-Rs, 1.5);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_sh_sl_db.AssertLessOrEqualsThan(-Rs, 1.5);    // Коэффициент передачи у верхнего края полосы подавления
        k_sh_db.AssertLessOrEqualsThan(-Rs, 0.36);      // Коэффициент передачи на верхней границе полосы подавления

        // Коэффициенты передачи в переходной полосе
        // между верхней полосой заграждения и верхней
        // полосой пропускания
        k_sh_ph_db.AssertLessThan(-Rp);
        k_ph_sh_db.AssertLessThan(-Rs);
        k_ph_sh_db.AssertLessThan(-Rp);

        // Коэффициенты передачи в верхней полосе пропускания
        // от верхней границы полосы пропускания до частоты,
        // близкой к половине частоты дискретизации
        k_ph_db.AssertLessThan(-Rp, 0.01);
        k_ph_fd05_db.AssertGreaterOrEqualsThan(-Rs);
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
        h_sl_pl.AssertLessOrEqualsThan(-Rp);

        h_sl.AssertLessOrEqualsThan(-Rs);
        h_sl_sh.AssertLessOrEqualsThan(-Rs);
        h_c0.AssertLessOrEqualsThan(-Rs);
        h_sh_sl.AssertLessOrEqualsThan(-Rs);
        h_sh.AssertLessOrEqualsThan(-Rs);

        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_ph_sh.AssertLessOrEqualsThan(-Rp);

        h_ph.AssertLessOrEqualsThan(-Rp);
        h_ph_fd05.AssertLessOrEqualsThan(-Rp);
        h_fd05.AssertGreaterOrEqualsThan(-Rp);
    }

    [TestMethod]
    public void TypeII_Corrected_Even_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.IICorrected);

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

        // Расчёт коэффициентов передачи по мощности
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

        // Расчёт коэффициентов передачи по мощности в логарифмическим масштабе
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
        // должны быть ниже уровня подавления Rs
        k_sl_db.AssertLessOrEqualsThan(-Rs, 2.94);      // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_sl_sh_db.AssertLessOrEqualsThan(-Rs);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertLessOrEqualsThan(-Rs, 0.4);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_sh_sl_db.AssertLessOrEqualsThan(-Rs, 0.5);    // Коэффициент передачи у верхнего края полосы подавления
        k_sh_db.AssertLessOrEqualsThan(-Rs, 1.1);       // Коэффициент передачи на верхней границе полосы подавления

        // Коэффициенты передачи в переходной полосе
        // между верхней полосой заграждения и верхней
        // полосой пропускания
        k_sh_ph_db.AssertLessThan(-Rp);
        k_ph_sh_db.AssertGreaterOrEqualsThan(-Rs);
        k_ph_sh_db.AssertLessOrEqualsThan(-Rp, 0.18);

        // Коэффициенты передачи в верхней полосе пропускания
        // от верхней границы полосы пропускания до частоты,
        // близкой к половине частоты дискретизации
        k_ph_db.AssertGreaterOrEqualsThan(-Rp);
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
        h_sh.AssertLessOrEqualsThan(-Rs, 9.5e-4);

        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_ph_sh.AssertLessOrEqualsThan(-Rp, 0.176).GreaterOrEqualsThan(-Rs);

        h_ph.AssertGreaterOrEqualsThan(-Rp);
        h_ph_fd05.AssertGreaterOrEqualsThan(-Rp);
        h_fd05.AssertGreaterOrEqualsThan(-Rp);
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_SignalProcessing()
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
        const double fph = 14 / Consts.pi2; // верхняя частота границы полосы пропускания

        var filter = new DSP.Filters.ChebyshevBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs, ChebyshevType.IICorrected);

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

        // Расчёт коэффициентов передачи по мощности
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

        // Расчёт коэффициентов передачи по мощности в логарифмическим масштабе
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
        // должны быть ниже уровня подавления Rs
        k_sl_db.AssertLessOrEqualsThan(-Rs, 3.6);       // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_sl_sh_db.AssertLessOrEqualsThan(-Rs, 1.27);   // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertLessOrEqualsThan(-Rs);            // Коэффициент передачи на центральной частоте полосы заграждения
        k_sh_sl_db.AssertLessOrEqualsThan(-Rs);         // Коэффициент передачи у верхнего края полосы подавления
        k_sh_db.AssertLessOrEqualsThan(-Rs, 1.4);       // Коэффициент передачи на верхней границе полосы подавления

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
        h_sh.AssertLessOrEqualsThan(-Rs, 1.3e-3);

        h_sh_ph.AssertLessOrEqualsThan(-Rp);
        h_ph_sh.AssertLessOrEqualsThan(-Rp).GreaterOrEqualsThan(-Rs);

        h_ph.AssertGreaterOrEqualsThan(-Rp);
        h_ph_fd05.AssertGreaterOrEqualsThan(-Rp);
        h_fd05.AssertGreaterOrEqualsThan(-Rp);
    }
}
