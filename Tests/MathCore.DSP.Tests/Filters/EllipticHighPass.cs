using System;
using System.Collections;
using System.Diagnostics;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;
using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions.EllipticJacobi;
// ReSharper disable ArgumentsStyleLiteral
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters;

[TestClassHandler("FailResultHandler")]
public class EllipticHighPass
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

        //switch (exception.Data["Expected"])
        //{
        //    case IEnumerable<Complex> expected:
        //        result.ToDebugEnum(expected);
        //        break;
        //    case IEnumerable expected:
        //        result.ToDebugEnum(expected);
        //        break;
        //    case { } expected:
        //        result.ToDebug(expected);
        //        break;
        //}
    }

    [TestMethod]
    public void Creation_Even()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        const double ws = Consts.pi2 * fs;
        const double wp = Consts.pi2 * fp;

        //var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);
        var Fp = DigitalFilter.ToDigitalFrequency(fp, dt);

        //(Fp, Fs).ToDebug();

        //var Ws = Consts.pi2 * Fs;
        var Wp = Consts.pi2 * Fp;

        var eps_p = Sqrt(Pow(10, Rp / 10) - 1);
        var eps_s = Sqrt(Pow(10, Rs / 10) - 1);

        const double kW = ws / wp;
        var k_eps = eps_p / eps_s;

        kW.AssertEquals(0.5);
        k_eps.AssertEquals(0.005088725841749188);

        var Kw = FullEllipticIntegral(kW);
        var Tw = FullEllipticIntegralComplimentary(kW);
        var K_eps = FullEllipticIntegral(k_eps);
        var T_eps = FullEllipticIntegralComplimentary(k_eps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N);             // порядок фильтра

        double_N.AssertEquals(3.317815612051342011);
        N.AssertEquals(4);

        var (L, r) = N.GetDivMod(2);

        var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);
        var m = (1 - k_eps * k_eps).Sqrt();

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));

        var k_W = Sqrt(1 - kp * kp);

        var v0_complex = sn_inverse((0, 1 / eps_p), k_eps) / N;

        var zeros = new Complex[N - r];
        var poles = new Complex[N];

        if (r != 0) poles[0] = Complex.i * sn_uk(v0_complex, k_W);
        for (var i = 0; i < L; i++)
        {
            // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_W);

            (poles[2 * i + r], poles[2 * i + r + 1]) = Complex.Conjugate(-p_re, p_im);
            (zeros[2 * i + 0], zeros[2 * i + 1]) = Complex.Conjugate(0, 1 / (k_W * cd_uk(u[i], k_W)));
        }

        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ +1.6095504012250093,
            /*[ 1]*/ -1.6095504012250093,
            /*[ 2]*/ +3.5252874329956083,
            /*[ 3]*/ -3.5252874329956083
            );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.105281264621164300, 0.993710811208774358),
            /*[ 1]*/ (-0.105281264621164300, -0.993710811208774358),
            /*[ 2]*/ (-0.364290595873426937, 0.478602767640667170),
            /*[ 3]*/ (-0.364290595873426937, -0.478602767640667170)
            );

        var high_pass_zeros_enum = AnalogBasedFilter.TransformToHighPassW(zeros, Wp);
        if (N.IsOdd())
            high_pass_zeros_enum = high_pass_zeros_enum.AppendFirst(0);
        var high_pass_zeros = high_pass_zeros_enum;
        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        high_pass_zeros.ToRe().Sum().AssertEquals(0);
        //high_pass_zeros.ToIm().ToDebugEnum();
        high_pass_zeros.ToIm().AssertEquals(
            /*[ 0]*/ -2.5188404830863624,
            /*[ 1]*/ +2.5188404830863624,
            /*[ 2]*/ -1.1500340857960620,
            /*[ 3]*/ +1.1500340857960620
        );

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/ (-0.427453184969549038, -4.034572083820466837),
            /*[ 1]*/ (-0.427453184969549038, 4.034572083820466837),
            /*[ 2]*/ (-4.082467720621100860, -5.363521243825417173),
            /*[ 3]*/ (-4.082467720621100860, 5.363521243825417173)
        );

        var z_zeros = DigitalFilter.ToZArray(high_pass_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(high_pass_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (0.968772524381009581, -0.247951196819950981),
            /*[ 1]*/ (0.968772524381009581, +0.247951196819950981),
            /*[ 2]*/ (0.993408901120038323, -0.114624409160865803),
            /*[ 3]*/ (0.993408901120038323, +0.114624409160865803)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.884631277392249560, -0.372228523605499628),
            /*[ 1]*/ (0.884631277392249560, 0.372228523605499628),
            /*[ 2]*/ (0.582466078525638697, -0.352438567686221504),
            /*[ 3]*/ (0.582466078525638697, 0.352438567686221504)
        );

        var k_zeros = z_zeros.Multiply(z => 1 + z).Re;
        var k_poles = z_poles.Multiply(z => 1 + z).Re;

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * k_poles / k_zeros;
        else
            g_norm = k_poles / k_zeros;

        g_norm.AssertEquals(0.550698192855287982);

        //var zz0 = Complex.Exp(-Consts.pi2 * 0 / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fs / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fp / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fd / 2 / fd);
        //var P0 = z_zeros.Multiply(z => 1 - z * zz0);
        //var Pp = z_poles.Multiply(z => 1 - z * zz0);
        //var k0 = g_norm * P0 / Pp;
        //k0.Abs.ToDebug();

        var B = GetCoefficientsInverted(z_zeros).ToArray(z => z.Re * g_norm);
        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        B.AssertEquals(
            /*[ 0]*/ +0.5506981928552880,
            /*[ 1]*/ -2.1611395301552800,
            /*[ 2]*/ +3.2213360608527040,
            /*[ 3]*/ -2.1611395301552796,
            /*[ 4]*/ +0.5506981928552878
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ +1,
            /*[ 1]*/ -2.9341947118357763,
            /*[ 2]*/ +3.4456770916870845,
            /*[ 3]*/ -1.8930671997108566,
            /*[ 4]*/ +0.4269234451315536
        );

        //Debug.WriteLine("---------------------------------------------");

        var filter = new DSP.Filters.EllipticHighPass(dt, fs, fp, Gp, Gs);

        filter.B.AssertEquals(Accuracy.Eps(1e-14), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-14), A);

        //var h_0 = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, 0, dt);
        //var h_p = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fs, dt);
        //var h_s = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fp, dt);
        //var h_d = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fd / 2, dt);

        //filter.GetTransmissionCoefficient(0).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fs).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fp).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fd/2).Abs.ToDebug();

        //var h_0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0, dt);
        //var h_p = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fs, dt);
        //var h_s = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fp, dt);
        //var h_d = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fd / 2, dt);

        //h_0.Abs.ToDebug();
        //h_p.Abs.ToDebug();
        //h_s.Abs.ToDebug();
        //h_d.Abs.ToDebug();
    }

    [TestMethod]
    public void Creation_Odd()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        const double ws = Consts.pi2 * fs;
        const double wp = Consts.pi2 * fp;

        //var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);
        var Fp = DigitalFilter.ToDigitalFrequency(fp, dt);

        //(Fp, Fs).ToDebug();

        //var Ws = Consts.pi2 * Fs;
        var Wp = Consts.pi2 * Fp;

        var eps_p = Sqrt(Pow(10, Rp / 10) - 1);
        var eps_s = Sqrt(Pow(10, Rs / 10) - 1);

        const double kW = ws / wp;
        var k_eps = eps_p / eps_s;

        kW.AssertEquals(0.75);
        k_eps.AssertEquals(0.005088725841749188);

        var Kw = FullEllipticIntegral(kW);
        var Tw = FullEllipticIntegralComplimentary(kW);
        var K_eps = FullEllipticIntegral(k_eps);
        var T_eps = FullEllipticIntegralComplimentary(k_eps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N);             // порядок фильтра

        double_N.AssertEquals(4.494923608977731355);
        N.AssertEquals(5);

        var (L, r) = N.GetDivMod(2);

        var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);
        var m = (1 - k_eps * k_eps).Sqrt();

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));

        var k_W = Sqrt(1 - kp * kp);

        var v0_complex = sn_inverse((0, 1 / eps_p), k_eps) / N;

        var zeros = new Complex[N - r];
        var poles = new Complex[N];

        if (r != 0) poles[0] = Complex.i * sn_uk(v0_complex, k_W);
        for (var i = 0; i < L; i++)
        {
            // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_W);

            (poles[2 * i + r], poles[2 * i + r + 1]) = Complex.Conjugate(-p_re, p_im);
            (zeros[2 * i + 0], zeros[2 * i + 1]) = Complex.Conjugate(0, 1 / (k_W * cd_uk(u[i], k_W)));
        }

        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ 1.253807568979517,
            /*[ 1]*/ -1.253807568979517,
            /*[ 2]*/ 1.7642884409084543,
            /*[ 3]*/ -1.7642884409084543
            );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/  -0.385344340275627639,
            /*[ 1]*/ (-0.049920708887154532, 0.998198050578146590),
            /*[ 2]*/ (-0.049920708887154532, -0.998198050578146590),
            /*[ 3]*/ (-0.219106729346206675, 0.741033961150657894),
            /*[ 4]*/ (-0.219106729346206675, -0.741033961150657894)
            );

        var high_pass_zeros_enum = AnalogBasedFilter.TransformToHighPassW(zeros, Wp);
        if (N.IsOdd())
            high_pass_zeros_enum = high_pass_zeros_enum.AppendFirst(0);
        var high_pass_zeros = high_pass_zeros_enum;
        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        high_pass_zeros.ToRe().Sum().AssertEquals(0);
        //high_pass_zeros.ToIm().ToDebugEnum();
        high_pass_zeros.ToIm().AssertEquals(
            /*[ 0]*/ +0,
            /*[ 1]*/ -3.2335111148461113,
            /*[ 2]*/ +3.2335111148461113,
            /*[ 3]*/ -2.297923976697309,
            /*[ 4]*/ +2.297923976697309
        );

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/  -10.520981590837891417,
            /*[ 1]*/ (-0.202613185262830536, -4.051386509914491008),
            /*[ 2]*/ (-0.202613185262830536, 4.051386509914491008),
            /*[ 3]*/ (-1.487597566404980531, -5.031156827179604996),
            /*[ 4]*/ (-1.487597566404980531, 5.031156827179604996)
        );

        var z_zeros = DigitalFilter.ToZArray(high_pass_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(high_pass_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/  1.000000000000000000,
            /*[ 1]*/ (0.949053713583807967, -0.315114342315266638),
            /*[ 2]*/ (0.949053713583807967, +0.315114342315266638),
            /*[ 3]*/ (0.973941725821161897, -0.226798401018385692),
            /*[ 4]*/ (0.973941725821161897, +0.226798401018385692)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/  0.310573838555953274,
            /*[ 1]*/ (0.903396072712712384, -0.381702758019330013),
            /*[ 2]*/ (0.903396072712712384, 0.381702758019330013),
            /*[ 3]*/ (0.764788580501989834, -0.413211764967296669),
            /*[ 4]*/ (0.764788580501989834, 0.413211764967296669)
        );

        var k_zeros = z_zeros.Multiply(z => 1 + z).Re;
        var k_poles = z_poles.Multiply(z => 1 + z).Re;

        double g_norm;
        if (N.IsEven())
            g_norm = Gp * k_poles / k_zeros;
        else
            g_norm = k_poles / k_zeros;

        g_norm.AssertEquals(0.5271810591547792);

        var zz0 = Complex.Exp(-Consts.pi2 * 0 / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fs / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fp / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fd / 2 / fd);
        //var P0 = z_zeros.Multiply(z => 1 - z * zz0);
        //var Pp = z_poles.Multiply(z => 1 - z * zz0);
        //var k0 = P0 / Pp;
        //k0.ToDebug();

        var B = GetCoefficientsInverted(z_zeros).ToArray(z => z.Re * g_norm);
        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        B.AssertEquals(
            /*[ 0]*/ 0.5271810591547792,
            /*[ 1]*/ -2.554714604145423,
            /*[ 2]*/ 5.031038000546882,
            /*[ 3]*/ -5.031038000546882,
            /*[ 4]*/ 2.5547146041454223,
            /*[ 5]*/ -0.5271810591547789
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ 1,
            /*[ 1]*/ -3.6469431449853578,
            /*[ 2]*/ 5.517284017908107,
            /*[ 3]*/ -4.228185429786587,
            /*[ 4]*/ 1.6077308828679036,
            /*[ 5]*/ -0.2257238521462117
        );

        var filter = new DSP.Filters.EllipticHighPass(dt, fs, fp, Gp, Gs);

        filter.B.AssertEquals(Accuracy.Eps(1e-14), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-14), A);
    }

    [TestMethod]
    public void TransmissionCoefficient_Even()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.EllipticHighPass(dt, fs, fp, Gp, Gs);

        var transmission__0 = filter.FrequencyResponse(0);
        var transmission_fs = filter.FrequencyResponse(fs);
        var transmission_fp = filter.FrequencyResponse(fp);
        var transmission_fd05 = filter.FrequencyResponse(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs, 1e-11);
        transmission_fs.Power.In_dB_byPower().AssertEquals(-Rs, 5e-3);
        transmission_fp.Power.In_dB_byPower().AssertEquals(-Rp, 6e-13);
        transmission_fd05.Power.In_dB_byPower().AssertEquals(-Rp, 1e-14);
    }

   [TestMethod]
    public void TransmissionCoefficient_Odd()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.EllipticHighPass(dt, fs, fp, Gp, Gs);

        var transmission__0 = filter.FrequencyResponse(0);
        var transmission_fs = filter.FrequencyResponse(fs);
        var transmission_fp = filter.FrequencyResponse(fp);
        var transmission_fd05 = filter.FrequencyResponse(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fs.Power.In_dB_byPower().AssertEquals(-Rs, 6.81e-1);
        transmission_fp.Power.In_dB_byPower().AssertEquals(-Rp, 8.9e-12);
        transmission_fd05.Power.In_dB_byPower().AssertEquals(0, 2e-15);
    }

    [TestMethod]
    public void SignalProcessing_Even()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.EllipticHighPass(dt, fs, fp, Gp, Gs);

        static EnumerableSignal GetSinSignal(double f0) => MathEnumerableSignal.Sin(dt, f0, (int)(100 * fd / (fs * 0.1)));

        var x_low = GetSinSignal(fs * 0.1);
        var x_s99 = GetSinSignal(fs * 0.99);
        var x_s = GetSinSignal(fs);

        var x_sp = GetSinSignal(fs + (fp - fs) * 0.1);
        var x_ps = GetSinSignal(fp - (fp - fs) * 0.1);

        var x_p = GetSinSignal(fp);

        var x_pd = GetSinSignal(fp + (fd / 2 - fp) * 0.1);
        var x_dp = GetSinSignal(fp + (fd / 2 - fp) * 0.9);

        var x_fd05 = GetSinSignal(0.9 * (fd / 2));

        /* ----------------------------------------------- */

        var y_low = filter.ProcessIndividual(x_low);
        var y_s99 = filter.ProcessIndividual(x_s99);
        var y_s = filter.ProcessIndividual(x_s);

        var y_sp = filter.ProcessIndividual(x_sp);
        var y_ps = filter.ProcessIndividual(x_ps);

        var y_p = filter.ProcessIndividual(x_p);

        var y_pd = filter.ProcessIndividual(x_pd);
        var y_dp = filter.ProcessIndividual(x_dp);

        var y_fd05 = filter.ProcessIndividual(x_fd05);

        /* ----------------------------------------------- */

        var k_low = y_low.Power / x_low.Power;
        var k_s99 = y_s99.Power / x_s99.Power;
        var k_s = y_s.Power / x_s.Power;

        var k_sp = y_sp.Power / x_sp.Power;
        var k_ps = y_ps.Power / x_ps.Power;

        var k_p = y_p.Power / x_p.Power;

        var k_pd = y_pd.Power / x_pd.Power;
        var k_dp = y_dp.Power / x_dp.Power;

        var k_fd05 = y_fd05.Power / x_fd05.Power;

        /* ----------------------------------------------- */

        var k_low_db = k_low.In_dB_byPower();
        var k_s99_db = k_s99.In_dB_byPower();

        var k_s_db = k_s.In_dB_byPower();

        var k_sp_db = k_sp.In_dB_byPower();
        var k_ps_db = k_ps.In_dB_byPower();

        var k_p_db = k_p.In_dB_byPower();

        var k_pd_db = k_pd.In_dB_byPower();
        var k_dp_db = k_dp.In_dB_byPower();

        var k_fd05_db = k_fd05.In_dB_byPower();

        /* --------------------------------------*/ //   -Rs  -Rp
                                                    // ---+----+->
        k_low_db.AssertLessThan(-Rs);               // \  |    |
        k_s99_db.AssertLessThan(-Rs, 2.7e-1);       //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 3e-1);   //   \+    | fs
                                                    //    \    |
        k_sp_db.AssertGreaterThan(-Rs, 0.5);        //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertGreaterThan(-Rs);             //       \ |
        k_ps_db.AssertLessThan(-Rp);                //        \+ fp
                                                    //         \0дБ
        k_p_db.AssertEquals(-Rp, 5.2e-3);           //          |
        k_pd_db.AssertGreaterOrEqualsThan(-Rp);     //          |
        k_dp_db.AssertGreaterOrEqualsThan(-Rp);     //          |
        k_fd05_db.AssertGreaterOrEqualsThan(-Rp);   //          |fd

        /* ----------------------------------------------- */

        var x =
            x_low +
            x_s99 +
            x_s +
            x_sp +
            x_ps +
            x_p +
            x_pd +
            x_dp +
            x_fd05;

        var y = filter.ProcessIndividual(x);

        var X = x.GetSpectrum();
        var Y = y.GetSpectrum();

        var H = Y / X;

        var h_low = H.GetValue(fs * 0.1).Power.In_dB_byPower();
        var h_s99 = H.GetValue(fs * 0.99).Power.In_dB_byPower();
        var h_s = H.GetValue(fs).Power.In_dB_byPower();

        var h_sp = H.GetValue(fs + (fp - fs) * 0.1).Power.In_dB_byPower();
        var h_ps = H.GetValue(fp - (fp - fs) * 0.1).Power.In_dB_byPower();

        var h_p = H.GetValue(fp).Power.In_dB_byPower();
        var h_pd = H.GetValue(fp + (fd / 2 - fp) * 0.1).Power.In_dB_byPower();
        var h_dp = H.GetValue(fp + (fd / 2 - fp) * 0.9).Power.In_dB_byPower();

        var h_fd05 = H.GetValue(0.9 * (fd / 2)).Power.In_dB_byPower();

        /* --------------------------------------*/ //   -Rs  -Rp
                                                    // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
                                                    //    \    |
        h_sp.AssertGreaterThan(-Rs, 1);             //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertGreaterThan(-Rs);                //       \ |
        h_ps.AssertLessThan(-Rp);                   //        \+ fp
                                                    //         \0дБ
        h_p.AssertEquals(-Rp, 8.1e-5);              //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }

    [TestMethod]
    public void SignalProcessing_Odd()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.EllipticHighPass(dt, fs, fp, Gp, Gs);

        static EnumerableSignal GetSinSignal(double f0) => MathEnumerableSignal.Sin(dt, f0, (int)(100 * fd / (fs * 0.1)));

        var x_low = GetSinSignal(fs * 0.1);
        var x_s99 = GetSinSignal(fs * 0.99);
        var x_s = GetSinSignal(fs);

        var x_sp = GetSinSignal(fs + (fp - fs) * 0.1);
        var x_ps = GetSinSignal(fp - (fp - fs) * 0.1);

        var x_p = GetSinSignal(fp);

        var x_pd = GetSinSignal(fp + (fd / 2 - fp) * 0.1);
        var x_dp = GetSinSignal(fp + (fd / 2 - fp) * 0.9);

        var x_fd05 = GetSinSignal(0.9 * (fd / 2));

        /* ----------------------------------------------- */

        var y_low = filter.ProcessIndividual(x_low);
        var y_s99 = filter.ProcessIndividual(x_s99);
        var y_s = filter.ProcessIndividual(x_s);

        var y_sp = filter.ProcessIndividual(x_sp);
        var y_ps = filter.ProcessIndividual(x_ps);

        var y_p = filter.ProcessIndividual(x_p);

        var y_pd = filter.ProcessIndividual(x_pd);
        var y_dp = filter.ProcessIndividual(x_dp);

        var y_fd05 = filter.ProcessIndividual(x_fd05);

        /* ----------------------------------------------- */

        var k_low = y_low.Power / x_low.Power;
        var k_s99 = y_s99.Power / x_s99.Power;
        var k_s = y_s.Power / x_s.Power;

        var k_sp = y_sp.Power / x_sp.Power;
        var k_ps = y_ps.Power / x_ps.Power;

        var k_p = y_p.Power / x_p.Power;

        var k_pd = y_pd.Power / x_pd.Power;
        var k_dp = y_dp.Power / x_dp.Power;

        var k_fd05 = y_fd05.Power / x_fd05.Power;

        /* ----------------------------------------------- */

        var k_low_db = k_low.In_dB_byPower();
        var k_s99_db = k_s99.In_dB_byPower();

        var k_s_db = k_s.In_dB_byPower();

        var k_sp_db = k_sp.In_dB_byPower();
        var k_ps_db = k_ps.In_dB_byPower();

        var k_p_db = k_p.In_dB_byPower();

        var k_pd_db = k_pd.In_dB_byPower();
        var k_dp_db = k_dp.In_dB_byPower();

        var k_fd05_db = k_fd05.In_dB_byPower();

        /* --------------------------------------*/ //   -Rs  -Rp
                                                    // ---+----+->
        k_low_db.AssertLessThan(-Rs);               // \  |    |
        k_s99_db.AssertLessThan(-Rs, 1.24);         //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 1.05);   //   \+    | fs
                                                    //    \    |
        k_sp_db.AssertGreaterThan(-Rs, 0.5);        //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertGreaterThan(-Rs);             //       \ |
        k_ps_db.AssertLessThan(-Rp);                //        \+ fp
                                                    //         \0дБ
        k_p_db.AssertEquals(-Rp, 1.6e-2);           //          |
        k_pd_db.AssertGreaterOrEqualsThan(-Rp);     //          |
        k_dp_db.AssertGreaterOrEqualsThan(-Rp);     //          |
        k_fd05_db.AssertGreaterOrEqualsThan(-Rp);   //          |fd

        /* ----------------------------------------------- */

        var x =
            x_low +
            x_s99 +
            x_s +
            x_sp +
            x_ps +
            x_p +
            x_pd +
            x_dp +
            x_fd05;

        var y = filter.ProcessIndividual(x);

        var X = x.GetSpectrum();
        var Y = y.GetSpectrum();

        var H = Y / X;

        var h_low = H.GetValue(fs * 0.1).Power.In_dB_byPower();
        var h_s99 = H.GetValue(fs * 0.99).Power.In_dB_byPower();
        var h_s = H.GetValue(fs).Power.In_dB_byPower();

        var h_sp = H.GetValue(fs + (fp - fs) * 0.1).Power.In_dB_byPower();
        var h_ps = H.GetValue(fp - (fp - fs) * 0.1).Power.In_dB_byPower();

        var h_p = H.GetValue(fp).Power.In_dB_byPower();
        var h_pd = H.GetValue(fp + (fd / 2 - fp) * 0.1).Power.In_dB_byPower();
        var h_dp = H.GetValue(fp + (fd / 2 - fp) * 0.9).Power.In_dB_byPower();

        var h_fd05 = H.GetValue(0.9 * (fd / 2)).Power.In_dB_byPower();

        /* --------------------------------------*/ //   -Rs  -Rp
                                                    // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
                                                    //    \    |
        h_sp.AssertLessOrEqualsThan(-Rs);           //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertGreaterThan(-Rs);                //       \ |
        h_ps.AssertLessThan(-Rp);                   //        \+ fp
                                                    //         \0дБ
        h_p.AssertEquals(-Rp, 8.1e-5);              //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }
}