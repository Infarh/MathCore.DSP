using System;
using System.Diagnostics;
using System.IO.Compression;
using System.Linq;
using System.Reflection;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Infrastructure;
using MathCore.DSP.Tests.Service;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Linq.Enumerable;
using static System.Math;

using static MathCore.SpecialFunctions;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ButterworthLowPass : UnitTest
{
    [TestMethod]
    public void Creation()
    {
        const double fd = 10;    // Гц // Частота дискретизации
        const double dt = 1 / fd; // 2с // Период дискретизации

        const double fp = 2 / Consts.pi2;   // Гц // Граничная частота полосы запирания
        const double fs = 4 / Consts.pi2;   // Гц // Граничная частота полосы пропускания

        Assert.IsTrue(fp < fs);
        Assert.IsTrue(fp < fd / 2);

        //const double wp = Consts.pi2 * fp * dt; // 0.628318530717959 рад/с
        //const double ws = Consts.pi2 * fs * dt; // 1.884955592153876 рад/с

        const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40; // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var eps_p = (10d.Pow(Rp / 10) - 1).Sqrt();
        var eps_s = (10d.Pow(Rs / 10) - 1).Sqrt();

        eps_p.AssertEquals(0.5088471399095875);
        eps_s.AssertEquals(99.99499987499375);

        var Fp = DigitalFilter.ToDigitalFrequency(fp, dt);  // Частота пропускания аналогового прототипа
        var Fs = DigitalFilter.ToDigitalFrequency(fs, dt);  // Частота подавления аналогового прототипа

        Fp.AssertEquals(0.319375180518077229);
        Fs.AssertEquals(0.645246083310777152);

        var Wp = Consts.pi2 * Fp;
        var Ws = Consts.pi2 * Fs;

        var kEps = eps_s / eps_p;
        var kW = Fs / Fp;

        kEps.AssertEquals(196.51284645671973);
        kW.AssertEquals(2.020338844941193202);

        var N = (int)Math.Ceiling(Log(kEps) / Log(kW));
        N.AssertEquals(8);

        var (L, r) = N.GetDivMod(2);

        var alpha = eps_p.Pow(-1d / N);
        alpha.AssertEquals(1.088119473662736647);

        var poles = new Complex[N];
        if (r != 0) poles[0] = -alpha;
        for (var (i, th0) = (r, Consts.pi05 / N); i < poles.Length; i += 2) 
            (poles[i], poles[i + 1]) = Complex.ConjugateExp(alpha, th0 * (i - r + 1 + N));

        //poles.ToDebugEnum();
        poles.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/ (-0.212281578508883212, +1.067211563088522608),
            /*[ 1]*/ (-0.212281578508883212, -1.067211563088522608),
            /*[ 2]*/ (-0.604526789535973275, +0.904738276905205363),
            /*[ 3]*/ (-0.604526789535973275, -0.904738276905205363),
            /*[ 4]*/ (-0.904738276905205363, +0.604526789535973497),
            /*[ 5]*/ (-0.904738276905205363, -0.604526789535973497),
            /*[ 6]*/ (-1.067211563088522608, +0.212281578508883656),
            /*[ 7]*/ (-1.067211563088522608, -0.212281578508883656)
        );

        var lowpass_poles = AnalogBasedFilter.TransformToLowPassW(poles, Wp);

        lowpass_poles.ToDebugEnum();
        lowpass_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/ (-0.425984051389412477, +2.141566444565760730),
            /*[ 1]*/ (-0.425984051389412477, -2.141566444565760730),
            /*[ 2]*/ (-1.213099943899241140, +1.815532366728786817),
            /*[ 3]*/ (-1.213099943899241140, -1.815532366728786817),
            /*[ 4]*/ (-1.815532366728786817, +1.213099943899241584),
            /*[ 5]*/ (-1.815532366728786817, -1.213099943899241584),
            /*[ 6]*/ (-2.141566444565760730, +0.425984051389413365),
            /*[ 7]*/ (-2.141566444565760730, -0.425984051389413365)
        );

        var z_poles = DigitalFilter.ToZArray(lowpass_poles, dt);

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(AccuracyComplex.Eps(1e-14),
            /*[ 0]*/ (0.936997507674962593, +0.203084896923829750),
            /*[ 1]*/ (0.936997507674962593, -0.203084896923829750),
            /*[ 2]*/ (0.871915749459332590, +0.160208721965230144),
            /*[ 3]*/ (0.871915749459332590, -0.160208721965230144),
            /*[ 4]*/ (0.827903822320540717, +0.101644552469968369),
            /*[ 5]*/ (0.827903822320540717, -0.101644552469968369),
            /*[ 6]*/ (0.805888478706402567, +0.034743688638417147),
            /*[ 7]*/ (0.805888478706402567, -0.034743688638417147)
        );

        var (g_norm, g_norm_im) = z_poles.Multiply(z => (1 - z) / 2);

        //g_norm.ToDebug();
        g_norm.AssertEquals(1.154301582524271E-08, 1.7e-024);
        g_norm_im.AssertEquals(0, 4.14e-25);

        var B = Range(0, N + 1).ToArray(i => g_norm * BinomialCoefficient(N, i));

        //B.ToDebugEnum();
        B.AssertEquals(Accuracy.Eps(1e-14),
            /*[ 0]*/ 1.154301582524271E-08,
            /*[ 1]*/ 9.234412660194168E-08,
            /*[ 2]*/ 3.232044431067959E-07,
            /*[ 3]*/ 6.464088862135918E-07,
            /*[ 4]*/ 8.080111077669897E-07,
            /*[ 5]*/ 6.464088862135918E-07,
            /*[ 6]*/ 3.232044431067959E-07,
            /*[ 7]*/ 9.234412660194168E-08,
            /*[ 8]*/ 1.154301582524271E-08
        );

        var (A, A_im) = Polynom.Array.GetCoefficientsInverted(z_poles);

        //A_im.Average(z => z * z).ToDebug();
        A_im.Average(z => z * z).AssertEquals(0, 1e-30);

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ +01,
            /*[ 1]*/ -06.8854111163224770,
            /*[ 2]*/ +20.8098097866188850,
            /*[ 3]*/ -36.0506821245047750,
            /*[ 4]*/ +39.1479188730887200,
            /*[ 5]*/ -27.2826323765932770,
            /*[ 6]*/ +11.9147694781671660,
            /*[ 7]*/ -02.9808064127161464,
            /*[ 8]*/ +00.3270368472739568
        );

        var filter = new DSP.Filters.ButterworthLowPass(dt, fp, fs, Gp, Gs);

        //filter.B.Zip(B, (fa, a) => fa - a).Max().ToDebug();

        filter.A.AssertEquals(A);
        filter.B.AssertEquals(Accuracy.Eps(1e-20), B);
    }

    [TestMethod]
    public void Creation_fp500_ps1500_fd_5000()
    {
        const double fp = 500;    // Гц // Граничная частота полосы пропускания
        const double fs = 1500;    // Гц // Граничная частота полосы запирания
        const double fd = 5000;      // Гц // Частота дискретизации
        const double dt = 1 / fd; // 2с // Период дискретизации
        const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
        const double Rs = 30; // Неравномерность в полосе пропускания (дБ)
        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var eps_p = (10d.Pow(Rp / 10) - 1).Sqrt();
        var eps_s = (10d.Pow(Rs / 10) - 1).Sqrt();

        var Fp = DigitalFilter.ToDigitalFrequency(fp, dt);  // Частота пропускания аналогового прототипа
        var Fs = DigitalFilter.ToDigitalFrequency(fs, dt);  // Частота подавления аналогового прототипа

        var k_eps = eps_s / eps_p;
        var k_W = Fs / Fp;

        var double_n = Log(k_eps) / Log(k_W);

        var N = (int)double_n;
        if (double_n - N > 0) N += 1;

        var L = N / 2;
        var r = N % 2;

        var alpha = eps_p.Pow(-1d / N);

        var th0 = PI / N;

        var poles = new Complex[N];
        if (r != 0) poles[0] = -alpha;
        for (var i = r; i < poles.Length; i += 2)
        {
            var w = th0 * (i + 1 - r - 0.5);
            poles[i] = (-alpha * Sin(w), alpha * Cos(w));
            poles[i + 1] = poles[i].ComplexConjugate;
        }

        var translated_poles = poles.ToArray(p => p * Consts.pi2 * Fp);
        var z_poles = translated_poles.ToArray(p => DigitalFilter.ToZ(p, dt));
        var kz = DigitalFilter.GetNormalizeCoefficient(translated_poles, dt);
        var WpN = (Consts.pi2 * Fp).Pow(N);
        var k = WpN * kz / eps_p;
        var b = Range(0, N + 1).ToArray(i => k * BinomialCoefficient(N, i));
        var a_complex = Polynom.Array.GetCoefficientsInverted(z_poles);
        var a = a_complex.ToRe() ?? throw new AssertFailedException("Отсутствует ссылка на массив вещественных значений");

        var filter = new DSP.Filters.ButterworthLowPass(dt, fp, fs, Gp, Gs);

        Assert.That.Value(filter.Order).IsEqual(N);

        var tr_k = filter.FrequencyResponse(fs, dt).Abs.In_dB();

        Assert.That.Collection(filter.A).IsEqualTo(a, 3e-12);
        Assert.That.Collection(filter.B).IsEqualTo(b, 1e-16);
    }

    [TestMethod]
    public void TransmissionCoefficient()
    {
        const double fd = 10;    // Гц // Частота дискретизации
        const double dt = 1 / fd; // 2с // Период дискретизации

        const double fp = 2 / Consts.pi2;   // Гц // Граничная частота полосы запирания
        const double fs = 4 / Consts.pi2;   // Гц // Граничная частота полосы пропускания

        const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40; // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ButterworthLowPass(dt, fp, fs, Gp, Gs);

        var transmission_0 = filter.FrequencyResponse(0, dt);
        var transmission_fp = filter.FrequencyResponse(fp, dt);
        var transmission_fs = filter.FrequencyResponse(fs, dt);
        var transmission_fd05 = filter.FrequencyResponse(fd / 2, dt);

        Assert.That.Value(transmission_0.Abs).IsEqual(1, 5.4e-10);
        Assert.That.Value(transmission_fp.Abs).IsEqual(Gp, 2.2e-10);
        Assert.That.Value(transmission_fs.Abs).LessOrEqualsThan(Gs, 1e-25);
        Assert.That.Value(transmission_fd05.Abs).IsEqual(0, 7e-26);
    }

    [TestMethod]
    public void TransmissionCoefficient_fp500_ps600_fd_5000()
    {
        const double fp = 500;    // Гц // Граничная частота полосы пропускания
        const double fs = 1500;    // Гц // Граничная частота полосы запирания
        const double fd = 5000;      // Гц // Частота дискретизации
        const double dt = 1 / fd; // 2с // Период дискретизации
        const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
        const double Rs = 30; // Неравномерность в полосе пропускания (дБ)
        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ButterworthLowPass(dt, fp, fs, Gp, Gs);

        var transmission_0 = filter.FrequencyResponse(0);
        var transmission_fp = filter.FrequencyResponse(fp);
        var transmission_fs = filter.FrequencyResponse(fs);
        var transmission_fd05 = filter.FrequencyResponse(fd / 2);

        Assert.That.Value(transmission_0.Abs).IsEqual(1, 5.56e-16);
        Assert.That.Value(transmission_fp.Abs).IsEqual(Gp, 8.9e-16);
        Assert.That.Value(transmission_fs.Abs).LessOrEqualsThan(Gs);
        Assert.That.Value(transmission_fd05.Abs).IsEqual(0, 9.61e-19);
    }

    [TestMethod]
    public void ImpulseResponse()
    {
        const double fp = 0.1;    // Гц // Граничная частота полосы пропускания
        const double fs = 0.3;    // Гц // Граничная частота полосы запирания
        const double fd = 1;      // Гц // Частота дискретизации
        const double dt = 1 / fd; // 2с // Период дискретизации
        const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
        const double Rs = 30; // Неравномерность в полосе пропускания (дБ)
        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ButterworthLowPass(dt, fp, fs, Gp, Gs);

        double[] expected_impulse_response =
        {
            +0.030466713814017596,
            +0.13656962408945683,
            +0.2655529821380315,
            +0.30340505314308713,
            +0.2307220811985055,
            +0.11400006843107224,
            +0.016215357391036463,
            -0.03502715933548784,
            -0.04382613004741935,
            -0.029116175888075326,
            -0.00954641846221929,
            +0.004003333255633864,
            +0.008889991723286804,
            +0.007517462792239407,
            +0.0036947291372435315,
            +0.0002967473365659246
        };
        var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();

        Assert.That.Collection(impulse_response)
           .IsEqualTo(expected_impulse_response, 1.7e-16);
    }

    [TestMethod]
    public void SignalProcessing()
    {
        const double fp = 0.1;    // Гц // Граничная частота полосы пропускания
        const double fs = 0.3;    // Гц // Граничная частота полосы запирания
        const double fd = 1;      // Гц // Частота дискретизации
        const double dt = 1 / fd; // 2с // Период дискретизации
        const double Rp = 1;      // Неравномерность в полосе пропускания (дБ)
        const double Rs = 30;     // Неравномерность в полосе пропускания (дБ)
        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ButterworthLowPass(dt, fp, fs, Gp, Gs);

        const int samples_count = 1024;
        // Сигнал s0(t) = 1
        var s_0 = new SamplesDigitalSignal(dt, Repeat(1d, samples_count));

        // Гармонические сигналы разной частоты и амплитудой = √2
        const double a0 = Consts.sqrt_2;
        // Сигнал с частотой равной частоте пропускания (граничной частоте фильтра)
        var s_fp = MathSamplesSignal.Cos(a0, fp, 0, dt, samples_count);
        // Сигнал с частотой равной частоте заграждения
        var s_fs = MathSamplesSignal.Cos(a0, fs, 0, dt, samples_count);
        // Сигнал с частотой равной половине частоты дискретизации
        var s_fd05 = MathSamplesSignal.Cos(a0, fd / 2, 0, dt, samples_count);

        var y_0 = filter.ProcessIndividual(s_0);
        var y_fp = filter.ProcessIndividual(s_fp);
        var y_fs = filter.ProcessIndividual(s_fs);
        var y_fd05 = filter.ProcessIndividual(s_fd05);

        // Постоянный сигнал не должен измениться своей мощности
        Assert.That.Value(y_0.Power).IsEqual(s_0.Power, 2.81e-3);
        // На граничной частоте сигнал должен быть ослаблен на коэффициент Gp (на Gp^2 по мощности)
        Assert.That.Value(y_fp.Power).IsEqual(s_0.Power * Gp * Gp, 3.2e-3);
        // На частоте заграждения сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
        Assert.That.Value(y_fs.Power).IsEqual(s_fs.Power * Gs * Gs, 1.84e-4);
        // На частоте в половину частоты дискретизации сигнал должен быть подавлен
        Assert.That.Value(y_fd05.Power).IsEqual(0, 1.33e-4);
    }
}