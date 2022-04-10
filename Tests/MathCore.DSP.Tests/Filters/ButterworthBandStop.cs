using MathCore.DSP.Filters;
using System;
using System.Diagnostics;
using System.Globalization;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;

using static MathCore.Polynom.Array;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ButterworthBandStop
{
    [TestMethod]
    public void Creation()
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

        var N = (int)Ceiling(Log(kEps) / Log(kW));
        N.AssertEquals(11); // порядок фильтра - число полюсов

        //var L = N / 2;  // число парных полюсов
        var r = N % 2;  // число (1 или 0) непарных (чисто действительных) полюсов

        var alpha = W0 * eps_p.Pow(-1d / N);
        alpha.AssertEquals(0.65401981150592248);

        var poles = new Complex[N];
        if (r != 0) poles[0] = -alpha;
        for (var (i, th0) = (r, Consts.pi05 / N); i < N; i += 2)
        {
            var w = th0 * (i - r + 1);
            var sin = -alpha * Sin(w);
            var cos = alpha * Cos(w);
            (poles[i], poles[i + 1]) = Complex.Conjugate(sin, cos);

            //var sin = -alpha * Sin(w);
            //var cos = alpha * Cos(w);
            
            //poles[i] = new Complex(sin, cos);
            //poles[i + 1] = new Complex(sin, -cos);
        }

        poles.AssertEquals(
            -0.6540198115059225,
            (-0.09307672370198979, +0.647362832843488),
            (-0.09307672370198979, -0.647362832843488),
            (-0.2716896485002241, +0.5949173461415183),
            (-0.2716896485002241, -0.5949173461415183),
            (-0.4282918937775254, +0.49427519416485316),
            (-0.4282918937775254, -0.49427519416485316),
            (-0.5501964769709404, +0.3535898055247178),
            (-0.5501964769709404, -0.3535898055247178),
            (-0.6275274137446106, +0.1842586737205135),
            (-0.6275274137446106, -0.1842586737205135)
        );

        var Hf = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(Enumerable.Empty<Complex>(), poles, F0);

        //var HHf = (W0.Power(N) / eps_p * Hf);
        //var HHf_abs = HHf.Abs;
        //var HHf_db = HHf.Power.In_dB_byPower();

        var sqrtWc = Wc.Sqrt();
        var pzf_zeros = Enumerable.Range(0, 2 * N).Select(i => i % 2 == 0 ? Complex.ImValue(sqrtWc) : Complex.ImValue(-sqrtWc));
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh).ToArray();

        pzf_poles.AssertEquals(
            (-7.3610426528900392, 1.1347289310788362),
            (-7.3610426528900392, -1.1347289310788362),
            (-0.3131230275290102, 3.1062866685921069),
            (-1.7820481618085897, -17.6785223734534398),
            (-0.3131230275290102, -3.1062866685921069),
            (-1.7820481618085897, 17.6785223734534398),
            (-0.9610736768408508, 3.0690350377114726),
            (-5.1547015818746615, -16.4607148701876405),
            (-0.9610736768408508, -3.0690350377114726),
            (-5.1547015818746615, 16.4607148701876405),
            (-1.6824629486202425, 2.9827065539257855),
            (-7.9584526399279962, -14.1089162573761602),
            (-1.6824629486202425, -2.9827065539257855),
            (-7.9584526399279962, 14.1089162573761602),
            (-2.5640115781587447, 2.8121729738581482),
            (-9.8209946959707644, -10.7715332082261313),
            (-2.5640115781587447, -2.8121729738581482),
            (-9.8209946959707644, 10.7715332082261313),
            (-3.8202228709510959, 2.4432367616673942),
            (-10.3055145368981247, -6.5909274969024576),
            (-3.8202228709510959, -2.4432367616673942),
            (-10.3055145368981247, 6.5909274969024576)
        );

        var h_F00 = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            0)
           .Power
           .In_dB_byPower();

        var h_Fpl = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fpl)
           .Power
           .In_dB_byPower();

        var h_Fsl = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fsl)
           .Power
           .In_dB_byPower();

        var Fc = (Fsl * Fsh).Sqrt();
        var h_Fc = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fc)
           .Power
           .In_dB_byPower();

        var h_Fsh = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fsh)
           .Power
           .In_dB_byPower();

        var h_Fph = DoubleArrayDSPExtensions.GetAnalogTransmissionCoefficientFromPoles(
            pzf_zeros,
            pzf_poles,
            Fph)
           .Power
           .In_dB_byPower();

        h_F00.AssertGreaterOrEqualsThan(-Rp);
        h_Fpl.AssertGreaterOrEqualsThan(-Rp);
        h_Fsl.AssertLessOrEqualsThan(-Rs);
        h_Fc.AssertLessOrEqualsThan(-Rs);
        h_Fsh.AssertLessOrEqualsThan(-Rs);
        h_Fph.AssertGreaterOrEqualsThan(-Rp);

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        z_zeros.AssertEquals(
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115),
            (0.75641755962253121, +0.65408904248175115),
            (0.75641755962253121, -0.65408904248175115)
        );

        z_poles.AssertEquals(
            (0.459422439986596431, 0.060525795245728167),
            (0.459422439986596431, -0.060525795245728167),
            (0.924174351815582495, 0.294245110857219050),
            (0.107109860752315725, -0.898541142586273822),
            (0.924174351815582495, -0.294245110857219050),
            (0.107109860752315725, 0.898541142586273822),
            (0.868248387868714389, 0.273541319968358398),
            (0.113391959227998079, -0.728580600328704753),
            (0.868248387868714389, -0.273541319968358398),
            (0.113391959227998079, 0.728580600328704753),
            (0.810546802219199058, 0.249064408686662936),
            (0.140304251667640267, -0.575441616956029023),
            (0.810546802219199058, -0.249064408686662936),
            (0.140304251667640267, 0.575441616956029023),
            (0.745619950387675923, 0.217558177990812951),
            (0.186530232602943796, -0.428582276793528394),
            (0.745619950387675923, -0.217558177990812951),
            (0.186530232602943796, 0.428582276793528394),
            (0.661762714414499520, 0.170446757573246371),
            (0.260282020315733054, -0.274089635120279951),
            (0.661762714414499520, -0.170446757573246371),
            (0.260282020315733054, 0.274089635120279951)
        );
           
        var G_norm = (N.IsOdd() ? 1 : Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        G_norm.AssertEquals(0.012857326566641241);

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
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
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fcc.Abs.AssertLessOrEqualsThan(Gs);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fd5.Abs.AssertGreaterOrEqualsThan(Gp);

        var filter = new DSP.Filters.ButterworthBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        filter.A.AssertEquals(A);
        filter.B.AssertEquals(B);
    }
}
