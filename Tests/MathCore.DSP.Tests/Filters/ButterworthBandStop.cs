using MathCore.DSP.Filters;

using static System.Math;

using static MathCore.Polynom.Array;
using MathCore.DSP.Signals;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ButterworthBandStop
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

        var Hf = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(Enumerable.Empty<Complex>(), poles, F0);

        //var HHf = (W0.Pow(N) / eps_p * Hf);
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

        var h_F00 = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            0)
           .Power
           .In_dB_byPower();

        var h_Fpl = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            Fpl)
           .Power
           .In_dB_byPower();

        var h_Fsl = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            Fsl)
           .Power
           .In_dB_byPower();

        var Fc = (Fsl * Fsh).Sqrt();
        var h_Fc = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            Fc)
           .Power
           .In_dB_byPower();

        var h_Fsh = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            Fsh)
           .Power
           .In_dB_byPower();

        var h_Fph = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
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
        h_fsl.Abs.AssertLessOrEqualsThan(Gs);
        h_fcc.Abs.AssertLessOrEqualsThan(Gs);
        h_fsh.Abs.AssertLessOrEqualsThan(Gs);
        h_fph.Abs.AssertGreaterOrEqualsThan(Gp);
        h_fd5.Abs.AssertGreaterOrEqualsThan(Gp);

        var filter = new DSP.Filters.ButterworthBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        filter.A.AssertEquals(Accuracy.Eps(1e-12), A);
        filter.B.AssertEquals(Accuracy.Eps(1e-12), B);
    }

    [TestMethod]
    public void TransmissionCoefficient()
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

        var filter = new DSP.Filters.ButterworthBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        var h_f0 = filter.FrequencyResponse(0);
        var h_pl = filter.FrequencyResponse(fpl);
        var h_sl = filter.FrequencyResponse(fsl);
        var h_c0 = filter.FrequencyResponse((fpl * fph).Sqrt());
        var h_sh = filter.FrequencyResponse(fsh);
        var h_ph = filter.FrequencyResponse(fph);
        var h_fd = filter.FrequencyResponse(fd / 2);

        var h_f0_db = h_f0.Power.In_dB_byPower();
        var h_pl_db = h_pl.Power.In_dB_byPower();
        var h_sl_db = h_sl.Power.In_dB_byPower();
        var h_c0_db = h_c0.Power.In_dB_byPower();
        var h_sh_db = h_sh.Power.In_dB_byPower();
        var h_ph_db = h_ph.Power.In_dB_byPower();
        var h_fd_db = h_fd.Power.In_dB_byPower();

        h_f0_db.AssertGreaterOrEqualsThan(-Rp);
        h_pl_db.AssertGreaterOrEqualsThan(-Rp);

        h_sl_db.AssertLessOrEqualsThan(-Rs);
        h_c0_db.AssertLessOrEqualsThan(-Rs);
        h_sh_db.AssertLessOrEqualsThan(-Rs);

        h_ph_db.AssertGreaterOrEqualsThan(-Rp);
        h_fd_db.AssertGreaterOrEqualsThan(-Rp);
    }

    [TestMethod]
    public void ImpulseResponse()
    {
        double[] expected_h =
        [
            +0.012857326529209, -0.078254837176753, +0.235536716014149, -0.405294467801370, +0.389957586557365,
            -0.094916570445630, -0.195017788724925, +0.166876192948785, +0.101455000886883, -0.136517871263717,
            -0.039418143498773, +0.138086585873978, +0.059945245448144, -0.071401263517821, -0.002748962419780,
            +0.118555916969853, +0.084923885237379, -0.005761180230300, +0.022730168856647, +0.113586876483959,
            +0.115515398529669, +0.050692136519340, +0.045532333951785, +0.103249475052354, +0.120356292899224,
            +0.074432968050876, +0.045804435029723, +0.068542913112179, +0.083635221028918, +0.051534860557429,
            +0.013199335297921, +0.010357111946136, +0.019045723066893, +0.001091969634829, -0.030746805041329,
            -0.041378115128919, -0.032860196189994, -0.033615982087794, -0.047464599803818, -0.051679901764387,
            -0.038361892405172, -0.024840246524572, -0.021973789316721, -0.018811221492894, -0.005212475303072,
            +0.010835204036393, +0.018134152406326, +0.019740862568349, +0.024568666420773, +0.031354624610029,
            +0.032151161571921, +0.026226069388840, +0.020175205325092, +0.016533679549206, +0.011191973002485,
            +0.002098748235656, -0.006840418859060, -0.012139580350658, -0.015372167233740, -0.019013829179943,
            -0.021806563938324, -0.021203027271530, -0.017696513500475, -0.013734139698559, -0.009904540335413,
            -0.004960103474912, +0.001119543743361, +0.006558333340733, +0.010238654906422, +0.012754977562679,
            +0.014679132451300, +0.015359926763873, +0.014102352937336, +0.011400347416904, +0.008214977642315,
            +0.004775256396496, +0.000866292552200, -0.003179230266005, -0.006557951365509, -0.008873758325497,
            -0.010304403336521, -0.010948406893043, -0.010565268912235, -0.009056823260525, -0.006771537897409,
            -0.004152037748306, -0.001377746959476, +0.001461393258891, +0.004075305558487, +0.006092281295104,
            +0.007354503620271, +0.007909261879345, +0.007780154470734, +0.006929646876150, +0.005432284237838,
            +0.003527001406808, +0.001463649835504, -0.000600428915935, -0.002526651937652, -0.004123667020692
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

        var filter = new DSP.Filters.ButterworthBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();

        Assert.That.Value(error2).LessThan(1.6e-14);
    }

    [TestMethod]
    public void SignalProcessing()
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

        var filter = new DSP.Filters.ButterworthBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

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