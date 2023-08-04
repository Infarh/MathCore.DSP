using System.Collections;
using System.Diagnostics;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions.EllipticJacobi;

namespace MathCore.DSP.Tests.Filters;

[TestClassHandler("FailResultHandler")]
public class EllipticBandStop : UnitTest
{
    // ReSharper disable once UnusedMember.Local
    private static void FailResultHandler(TestResult result)
    {
        if (result.TestFailureException?.InnerException is not AssertFailedException exception) return;
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
        const double W1 = 1;                        // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var F0 = W0 / Consts.pi2;
        //const double F1 = 1 / Consts.pi2;

        W0.AssertEquals(0.615059351152204);
        F0.AssertEquals(0.0978897360307671);

        var eps_p = Sqrt(10.Power(Rp / 10) - 1);
        var eps_s = Sqrt(10.Power(Rs / 10) - 1);

        eps_p.AssertEquals(0.50884713990958752);
        eps_s.AssertEquals(99.994999874993752);

        var kW = W0 / W1;
        var k_eps = eps_p / eps_s;

        kW.AssertEquals(0.615059351152204);

        var Kw = FullEllipticIntegral(kW);
        var Tw = FullEllipticIntegralComplimentary(kW);
        var K_eps = FullEllipticIntegral(k_eps);
        var T_eps = FullEllipticIntegralComplimentary(k_eps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N);             // порядок фильтра

        double_N.AssertEquals(3.7906792606389264);
        N.AssertEquals(4);

        var L = N / 2;  // число парных полюсов
        var r = N % 2;  // 1, или 0 - число непарных полюсов

        var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);
        var m = (1 - k_eps * k_eps).Sqrt();

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));

        var k_W = Sqrt(1 - kp * kp);

        var im_pz = Enumerable.Range(0, L).ToArray(i => 1 / (k_W * cd_uk(u[i], k_W)));

        var v0_complex = sn_inverse((0, 1 / eps_p), k_eps) / N;

        var zeros = new Complex[N - r];
        var poles = new Complex[N];

        if (r != 0) poles[0] = Complex.ImValue(W0) * sn_uk(v0_complex, k_W);
        for (var i = 0; i < L; i++)
        {
            // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_W);

            var pp = new Complex(-p_re * W0, p_im * W0);

            poles[r + 2 * i] = pp;
            poles[r + 2 * i + 1] = pp.ComplexConjugate;

            var p0_im = W0 / (k_W * cd_uk(u[i], k_W));
            zeros[2 * i] = (0, p0_im);
            zeros[2 * i + 1] = zeros[2 * i].ComplexConjugate;
        }

        zeros.AssertEquals(
            (0, +0.98996902542422383),
            (0, -0.98996902542422383),
            (0, +2.1682610011632977),
            (0, -2.1682610011632977)
            );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.064754226306376811, 0.611191126774999094),
            /*[ 1]*/ (-0.064754226306376811, -0.611191126774999094),
            /*[ 2]*/ (-0.224060337528759729, 0.294369107724717805),
            /*[ 3]*/ (-0.224060337528759729, -0.294369107724717805)
            );

        var pzf_zeros = AnalogBasedFilter.TransformToBandStop(zeros, Fsl, Fsh);
        var pzf_poles = AnalogBasedFilter.TransformToBandStop(poles, Fsl, Fsh);

        var wc_sqrt = wc.Sqrt();
        if (zeros.Length != poles.Length)
            pzf_zeros = pzf_zeros.AppendLast(
                (0, +wc_sqrt),
                (0, -wc_sqrt));

        pzf_zeros.AssertEquals(
            (+5.44664340994985E-16, +4.0319948722574663),
            (-5.44664340994985E-16, -13.758092567621574),
            (+5.44664340994985E-16, +13.758092567621574),
            (-5.44664340994985E-16, -4.0319948722574663),

            (+4.758917038943161E-16, +5.551565428550191),
            (-4.758917038943161E-16, -9.9922372164459166),
            (+4.758917038943161E-16, +9.9922372164459166),
            (-4.758917038943161E-16, -5.551565428550191)
            );

        //pzf_poles.ToDebugEnum();
        pzf_poles.AssertEquals(
            /*[ 0]*/ (-0.227955413152198361, 2.972703475040113119),
            /*[ 1]*/ (-1.422586359987087601, -18.551555137033666654),
            /*[ 2]*/ (-0.227955413152198361, -2.972703475040113119),
            /*[ 3]*/ (-1.422586359987087601, 18.551555137033666654),
            /*[ 4]*/ (-1.130719053614130054, 1.734335688271796627),
            /*[ 5]*/ (-14.633073912059700206, -22.444710941843673879),
            /*[ 6]*/ (-1.130719053614130054, -1.734335688271796627),
            /*[ 7]*/ (-14.633073912059700206, 22.444710941843673879)
        );

        var h_F00 = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            0)
           .Abs;

        var h_Fpl = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            Fpl)
           .Abs;

        var h_Fsl = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            Fsl)
            .Abs;

        var Fc = (Fsl * Fsh).Sqrt();
        var h_Fc = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            Fc)
           .Abs;

        var h_Fsh = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            Fsh).Abs;

        var h_Fph = DoubleArrayDSPExtensions.AnalogFrequencyResponseFromPoles(
            pzf_zeros,
            pzf_poles,
            Fph)
           .Abs;

        h_F00.AssertGreaterOrEqualsThan(Gp);
        h_Fpl.AssertGreaterOrEqualsThan(Gp);
        h_Fsl.AssertLessOrEqualsThan(Gs);
        h_Fc.AssertLessOrEqualsThan(Gs, 1e-2);
        h_Fsh.AssertLessOrEqualsThan(Gs, 1e-3);
        h_Fph.AssertGreaterOrEqualsThan(Gp);

        var z_zeros = DigitalFilter.ToZArray(pzf_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(pzf_poles, dt);

        z_zeros.AssertEquals(
            (0.92188968196320531, +0.38745246713600884),
            (0.35757714717701156, -0.9338836029274471),
            (0.35757714717701156, +0.9338836029274471),
            (0.92188968196320531, -0.38745246713600884),

            (0.85692452818129894, +0.51544190070390872),
            (0.60049677950867342, -0.79962717425042007),
            (0.60049677950867342, +0.79962717425042007),
            (0.85692452818129894, -0.51544190070390872)
        );

        z_poles.AssertEquals(
            /*[ 0]*/ (0.935656421150197093, 0.284464368845482374),
            /*[ 1]*/ (0.067011448261103043, -0.924011759440699354),
            /*[ 2]*/ (0.935656421150197093, -0.284464368845482374),
            /*[ 3]*/ (0.067011448261103043, 0.924011759440699354),
            /*[ 4]*/ (0.880311827269883440, 0.154329433789710341),
            /*[ 5]*/ (-0.186642278225500613, -0.527114024123290781),
            /*[ 6]*/ (0.880311827269883440, -0.154329433789710341),
            /*[ 7]*/ (-0.186642278225500613, 0.527114024123290781)
            );

        var G_norm = (N.IsOdd() ? 1 : Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var h_f00 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0, dt).Abs.AssertGreaterOrEqualsThan(Gp);
        var h_fpl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fpl, dt).Abs.AssertGreaterOrEqualsThan(Gp);
        var h_fsl = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsl, dt).Abs.AssertLessOrEqualsThan(Gs);
        var h_fc = DoubleArrayDSPExtensions.FrequencyResponse(A, B, (fsl * fsh).Sqrt(), dt).Abs.AssertLessOrEqualsThan(Gs);
        var h_fsh = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fsh, dt).Abs.AssertLessOrEqualsThan(Gs);
        var h_fph = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fph, dt).Abs.AssertGreaterOrEqualsThan(Gp);
        var h_fd5 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fd / 2, dt).Abs.AssertGreaterOrEqualsThan(Gp);


        var filter = new DSP.Filters.EllipticBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        var filter_A = filter.A;
        var filter_B = filter.B;

        Assert.That.Collection(filter_A).IsEqualTo(A, 1.156e-14);
        Assert.That.Collection(filter_B).IsEqualTo(B, 1.156e-14);


        //var h_f0 = filter.GetTransmissionCoefficient(0).Power.In_dB_byPower();
        //var h_sl = filter.GetTransmissionCoefficient(fsl).Power.In_dB_byPower();
        //var h_pl = filter.GetTransmissionCoefficient(fpl).Power.In_dB_byPower();
        //var h_c0 = filter.GetTransmissionCoefficient((fpl * fph).Sqrt()).Power.In_dB_byPower();
        //var h_ph = filter.GetTransmissionCoefficient(fph).Power.In_dB_byPower();
        //var h_sh = filter.GetTransmissionCoefficient(fsh).Power.In_dB_byPower();
        //var h_fd = filter.GetTransmissionCoefficient(fd / 2).Power.In_dB_byPower();

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

        var filter = new DSP.Filters.EllipticBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

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

        h_ph_db.AssertGreaterOrEqualsThan(-Rp, 1e-14);
        h_fd_db.AssertGreaterOrEqualsThan(-Rp, 1e-14);
    }

    [TestMethod]
    public void ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.183172689114915, -0.381198732755215, +0.453076019965961, -0.062912057693868, -0.148022513286273, +0.181720917354020,
            +0.159684282091557, -0.008939212573129, +0.013823995278017, +0.135554059168909, +0.148284929774527, +0.056004332116796,
            +0.025548260646322, +0.094885335574295, +0.122424902141717, +0.052056841966505, +0.001065090591171, +0.034442621318458,
            +0.060123873245061, +0.011354378464955, -0.040356296485409, -0.027077829275585, +0.000044420317663, -0.023128237897633,
            -0.060594561288147, -0.052073449389965, -0.020287979095019, -0.020266904487746, -0.040437370172537, -0.032826559076297,
            -0.001936058603854, +0.009350356552149, -0.001974752761781, -0.000562321247940, +0.019999566957344, +0.029535095680068,
            +0.018056789899852, +0.009963027702425, +0.017209104658635, +0.020821399249043, +0.008568929474867, -0.004180055131952,
            -0.004063764912808, -0.001667220690710, -0.009121510302649, -0.018330967469683, -0.017328894005210, -0.011188193780463,
            -0.010799046197802, -0.013412032494572, -0.009639736372143, -0.001057609445802, +0.003494178337189, +0.003566818086069,
            +0.006082217413436, +0.011810433583898, +0.014477215040405, +0.012283432231163, +0.010234382169118, +0.010569981365131,
            +0.009473710757618, +0.004823124654922, +0.000020902414037, -0.002183356432381, -0.003852259835463, -0.007157959668392,
            -0.010220735160551, -0.010618965224709, -0.009453240708199, -0.008889556020284, -0.008341509961118, -0.006056033020880,
            -0.002581994165077, +0.000153890697998, +0.001980517842253, +0.004143387218973, +0.006602558130809, +0.007990170563372,
            +0.007946215986297, +0.007482399464124, +0.007031634927994, +0.005870874162197, +0.003706126172102, +0.001361646805830,
            -0.000511736057887, -0.002227698614659, -0.004088678238406, -0.005601108583327, -0.006228765259648, -0.006182802571738,
            -0.005889948947828, -0.005244962394558, -0.003961757092649, -0.002250752695381, -0.000592009948249, +0.000896174847138,
            +0.002368363895985, +0.003718994249984, +0.004606024681828, +0.004926740575332
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

        var filter = new DSP.Filters.EllipticBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2() / 2).Sum();

        Assert.That.Value(error2).LessThan(1e-10);
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

        var filter = new DSP.Filters.EllipticBandStop(dt, fpl, fsl, fsh, fph, Gp, Gs);

        //const double total_time = 1 / fpl;
        //const int samples_count = (int)(total_time * fd) + 1;

        // Метод формирования гармонического сигнала на заданной частоте
        static EnumerableSignal GetSinSignal(double f0) => MathEnumerableSignal.Sin(dt, f0, (int)(100 * fd / (fpl * 0.1)));

        var x_low     = GetSinSignal(fpl * 0.1);                    // частота, близкая к нулю
        var x_pl99    = GetSinSignal(fpl * 0.99);                   // частота чуть ниже нижней границы пропускания
        var x_pl      = GetSinSignal(fpl);                          // частота на нижней границе пропускания

        var x_pl_sl   = GetSinSignal(fpl + (fsl - fpl) * 0.1);      // частота выше нижней границы пропускания
        var x_sl_pl   = GetSinSignal(fsl - (fsl - fpl) * 0.1);      // частота ниже нижней границы подавления

        var x_sl      = GetSinSignal(fsl);                          // частота на границе подавления
        var x_sl_sh   = GetSinSignal(fsl + (fsh - fsl) * 0.1);      // частота выше нижней границы подавления
        var x_c0      = GetSinSignal((fsl * fsh).Sqrt());           // частота в середине полосы подавления
        var x_sh_sl   = GetSinSignal(fsh - (fsh - fsl) * 0.1);      // частота ниже верхней границы подавления
        var x_sh      = GetSinSignal(fsh);                          // частота на границе подавления

        var x_sh_ph   = GetSinSignal(fsh + (fph - fsh) * 0.1);      // частота выше верхней границы подавления
        var x_ph_sh   = GetSinSignal(fph - (fph - fsh) * 0.1);      // частота ниже верхней границы пропускания

        var x_ph      = GetSinSignal(fph);                          // частота на верхней границе пропускания
        var x_ph_fd05 = GetSinSignal(fph + (fd / 2 - fph) * 0.1);   // частота выше верхней границы пропускания
        var x_fd05    = GetSinSignal(0.9 * (fd / 2));               // частота ниже половины частоты дискретизации

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
        k_sl_db.AssertLessOrEqualsThan(-Rs);            // Коэффициент передачи на нижней границе полосы заграждения должен быть не больше Rs (-40дБ и ниже)
        k_sl_sh_db.AssertLessOrEqualsThan(-Rs);         // Коэффициент передачи у нижнего края полосы заграждения должен быть меньше Rs
        k_c0_db.AssertLessOrEqualsThan(-Rs, 1.5);       // Коэффициент передачи на центральной частоте полосы заграждения
        k_sh_sl_db.AssertLessOrEqualsThan(-Rs, 1.5);    // Коэффициент передачи у верхнего края полосы подавления
        k_sh_db.AssertLessOrEqualsThan(-Rs);            // Коэффициент передачи на верхней границе полосы подавления

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
        var h_low     = H.GetValue(fpl * 0.1).Power.In_dB_byPower();
        var h_pl99    = H.GetValue(fpl * 0.99).Power.In_dB_byPower();
        var h_pl      = H.GetValue(fpl).Power.In_dB_byPower();

        var h_pl_sl   = H.GetValue(fpl + (fsl - fpl) * 0.1).Power.In_dB_byPower();
        var h_sl_pl   = H.GetValue(fsl - (fsl - fpl) * 0.1).Power.In_dB_byPower();

        var h_sl      = H.GetValue(fsl).Power.In_dB_byPower();
        var h_sl_sh   = H.GetValue(fsl + (fsh - fsl) * 0.1).Power.In_dB_byPower();
        var h_c0      = H.GetValue((fsl * fsh).Sqrt()).Power.In_dB_byPower();
        var h_sh_sl   = H.GetValue(fsh - (fsh - fsl) * 0.1).Power.In_dB_byPower();
        var h_sh      = H.GetValue(fsh).Power.In_dB_byPower();
                      
        var h_sh_ph   = H.GetValue(fsh + (fph - fsh) * 0.1).Power.In_dB_byPower();
        var h_ph_sh   = H.GetValue(fph - (fph - fsh) * 0.1).Power.In_dB_byPower();
                      
        var h_ph      = H.GetValue(fph).Power.In_dB_byPower();
        var h_ph_fd05 = H.GetValue(fph + (fd / 2 - fph) * 0.1).Power.In_dB_byPower();
        var h_fd05    = H.GetValue(0.9 * (fd / 2)).Power.In_dB_byPower();

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