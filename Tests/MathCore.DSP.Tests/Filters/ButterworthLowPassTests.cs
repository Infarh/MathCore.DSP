using System;
using System.Linq;
using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class ButterworthLowPassTests : UnitTest
    {
        [TestMethod]
        public void CreationTest()
        {
            const double fd = 0.5;    // Гц // Частота дискретизации
            const double dt = 1 / fd; // 2с // Период дискретизации

            const double fp = 0.05;    // Гц // Граничная частота полосы пропускания
            const double fs = 0.15;    // Гц // Граничная частота полосы запирания

            Assert.IsTrue(fp < fs);
            Assert.IsTrue(fp < fd / 2);

            //const double wp = Consts.pi2 * fp * dt; // 0.628318530717959 рад/с
            //const double ws = Consts.pi2 * fs * dt; // 1.884955592153876 рад/с

            const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 30; // Неравномерность в полосе пропускания (дБ)

            var eps_p = (10d.Pow(Rp / 10) - 1).Sqrt();
            var eps_s = (10d.Pow(Rs / 10) - 1).Sqrt();

            Assert.That.Value(eps_p).IsEqual(0.508847139909588, 4.441e-16);
            Assert.That.Value(eps_s).IsEqual(31.606961258558215);

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            Assert.That.Value(Gp).IsEqual(0.891250938133746, 4.45e-16);
            Assert.That.Value(Gs).IsEqual(0.0316227766016838, 4.45e-16);

            Assert.That.Value(0.891250938133746.In_dB()).IsEqual(-1, 4.67e-15);
            Assert.That.Value(0.0316227766016838.In_dB()).IsEqual(-30);

            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);  // Частота подавления аналогового протипа

            Assert.That.Value(Fp).IsEqual(0.051712575763384, 2.5e-16);
            Assert.That.Value(Fs).IsEqual(0.219057986225303, 3.89e-16);

            var Wp = Consts.pi2 * Fp;
            var Ws = Consts.pi2 * Fs;

            var k_eps = eps_s / eps_p;
            var k_W = Ws / Wp;

            Assert.That.Value(k_eps).GreaterThan(0);
            Assert.That.Value(k_W).GreaterThan(0);

            var double_n = Math.Log(k_eps) / Math.Log(k_W);
            // Порядок фильтра
            var N = (int)double_n;
            if (double_n - N > 0) N += 1;

            Assert.That.Value(N).IsEqual(3);

            var L = N / 2;
            var r = N % 2;

            Assert.That.Value(L).IsEqual(1);
            Assert.That.Value(r).IsEqual(1);
            Assert.That.Value(2 * L + r).IsEqual(N);

            var alpha = Math.Pow(eps_p, -1d / N);
            Assert.That.Value(alpha).IsEqual(1.252576388181026, 4.45e-16);

            var th0 = Math.PI / N;

            var poles = new Complex[N];
            if (r != 0) poles[0] = -alpha;
            for (var i = r; i < poles.Length; i += 2)
            {
                var w = th0 * (i + 1 - r - 0.5);
                var sin = -alpha * Math.Sin(w);
                var cos = alpha * Math.Cos(w);
                poles[i] = new Complex(sin, cos);
                poles[i + 1] = new Complex(sin, -cos);
            }

            Assert.That.Value(poles[0].Re).IsEqual(-1.252576388181026, 4.45e-16);
            Assert.That.Value(poles[0].Im).IsEqual(0);

            Assert.That.Value(poles[1].Re).IsEqual(-0.626288194090513, 4.45e-16);
            Assert.That.Value(poles[1].Im).IsEqual(1.084762972345327, 4.45e-16);

            Assert.That.Value(poles[2].Re).IsEqual(-0.626288194090513, 4.45e-16);
            Assert.That.Value(poles[2].Im).IsEqual(-1.084762972345327, 4.45e-16);

            var translated_poles = poles.ToArray(p => p * Wp);

            Assert.That.Value(translated_poles[0].Re).IsEqual(-0.40698673955629);
            Assert.That.Value(translated_poles[0].Im).IsEqual(0);

            Assert.That.Value(translated_poles[1].Re).IsEqual(-0.203493369778145, 1.12e-16);
            Assert.That.Value(translated_poles[1].Im).IsEqual(0.352460855459148, 4.45e-16);

            Assert.That.Value(translated_poles[2].Re).IsEqual(-0.203493369778145, 1.12e-16);
            Assert.That.Value(translated_poles[2].Im).IsEqual(-0.352460855459148, 4.45e-16);

            var z_poles = translated_poles.ToArray(p => DigitalFilter.ToZ(p, dt));

            Assert.That.Value(z_poles[0].Re).IsEqual(0.421477504919999, 2.23e-16);
            Assert.That.Value(z_poles[0].Im).IsEqual(0);

            Assert.That.Value(z_poles[1].Re).IsEqual(0.530553579281761, 2.23e-16);
            Assert.That.Value(z_poles[1].Im).IsEqual(0.4482452811345, 4.45e-16);

            Assert.That.Value(z_poles[2].Re).IsEqual(0.530553579281761, 2.23e-16);
            Assert.That.Value(z_poles[2].Im).IsEqual(-0.4482452811345, 4.45e-16);

            Assert.That.Value(z_poles[0].Abs).LessThan(1);
            Assert.That.Value(z_poles[1].Abs).LessThan(1);
            Assert.That.Value(z_poles[2].Abs).LessThan(1);

            var kz = DigitalFilter.GetNormalizeCoefficient(translated_poles, dt);
            Assert.That.Value(kz).IsEqual(0.451944218734017, 1.12e-16);

            var WpN = Wp.Pow(N);
            Assert.That.Value(WpN).IsEqual(0.034302685030762, 2.78e-16);

            var k = WpN * kz / eps_p;
            Assert.That.Value(k).IsEqual(0.030466713814018, 4.1e-16);

            var b = new double[N + 1];
            for (var i = 0; i < b.Length; i++)
                b[i] = k * SpecialFunctions.BinomialCoefficient(N, i);

            Assert.That.Value(b[0]).IsEqual(1 * k);
            Assert.That.Value(b[1]).IsEqual(3 * k);
            Assert.That.Value(b[2]).IsEqual(3 * k);
            Assert.That.Value(b[3]).IsEqual(1 * k);

            var a_complex = Polynom.Array.GetCoefficientsInverted(z_poles);
            Assert.That.Value(a_complex.Length).IsEqual(4);

            Assert.That.Value(a_complex[0].Im).IsEqual(0);
            Assert.That.Value(a_complex[1].Im).IsEqual(0);
            Assert.That.Value(a_complex[2].Im).IsEqual(0, 2.78e-1);
            Assert.That.Value(a_complex[3].Im).IsEqual(0);

            var a = a_complex.ToRe() ?? throw new AssertFailedException("Отсутствует ссылка на массив вещественных значений");

            Assert.That.Value(a[0]).IsEqual(1);
            Assert.That.Value(a[1]).IsEqual(-1.482584663483521, 6.67e-16);
            Assert.That.Value(a[2]).IsEqual(0.929643730192138, 1.12e-16);
            Assert.That.Value(a[3]).IsEqual(-0.203325356196476, 3.1e-16);

            var filter = new ButterworthLowPass(fp, fs, dt, Gp, Gs);

            var A = filter.A;
            var B = filter.B;

            Assert.That.Value(A.Count).IsEqual(a.Length);
            Assert.That.Value(B.Count).IsEqual(b.Length);

            Assert.That.Value(A[0]).IsEqual(a[0], 5.56e-17);
            Assert.That.Value(A[1]).IsEqual(a[1], 2.23e-16);
            Assert.That.Value(A[2]).IsEqual(a[2], 5.56e-17);
            Assert.That.Value(A[3]).IsEqual(a[3], 5.56e-17);

            Assert.That.Value(B[0]).IsEqual(b[0], 1.39e-17);
            Assert.That.Value(B[1]).IsEqual(b[1], 2.78e-17);
            Assert.That.Value(B[2]).IsEqual(b[2], 2.78e-17);
            Assert.That.Value(B[3]).IsEqual(b[3], 1.39e-17);
        }

        [TestMethod]
        public void TransmissionCoefficientTest()
        {
            const double fp = 0.1;    // Гц // Граничная частота полосы пропускания
            const double fs = 0.3;    // Гц // Граничная частота полосы запирания
            const double fd = 1;      // Гц // Частота дискретизации
            const double dt = 1 / fd; // 2с // Период дискретизации
            const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 30; // Неравномерность в полосе пропускания (дБ)
            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ButterworthLowPass(fp, fs, dt, Gp, Gs);

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            Assert.That.Value(transmission_0.Abs).IsEqual(1, 2.23e-16);
            Assert.That.Value(transmission_fp.Abs).GreaterOrEqualsThan(Gp);
            Assert.That.Value(transmission_fs.Abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05.Abs).IsEqual(0, 4.27e-34);
        }

        [TestMethod]
        public void ImpulseResponseTest()
        {
            const double fp = 0.1;    // Гц // Граничная частота полосы пропускания
            const double fs = 0.3;    // Гц // Граничная частота полосы запирания
            const double fd = 1;      // Гц // Частота дискретизации
            const double dt = 1 / fd; // 2с // Период дискретизации
            const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 30; // Неравномерность в полосе пропускания (дБ)
            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ButterworthLowPass(fp, fs, dt, Gp, Gs);

            double[] expected_impulse_response =
            {
                +0.030466713814017603,
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
               .IsEqualTo(expected_impulse_response);
        }

        [TestMethod]
        public void SignallProcessingTest()
        {
            const double fp = 0.1;    // Гц // Граничная частота полосы пропускания
            const double fs = 0.3;    // Гц // Граничная частота полосы запирания
            const double fd = 1;      // Гц // Частота дискретизации
            const double dt = 1 / fd; // 2с // Период дискретизации
            const double Rp = 1;      // Неравномерность в полосе пропускания (дБ)
            const double Rs = 30;     // Неравномерность в полосе пропускания (дБ)
            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ButterworthLowPass(fp, fs, dt, Gp, Gs);

            const int samples_count = 1024;
            // Сигнал s0(t) = 1
            var s_0 = new SamplesDigitalSignal(dt, Enumerable.Repeat(1d, samples_count));

            // Гармонические сигналы разной частоты и амплитудой = √2
            const double a0 = Consts.sqrt_2;
            // Сигнал с частотой равной частоте пропускания (граничной частоте фильтра)
            var s_fp = new SamplesDigitalSignal(dt, samples_count, t => a0 * Math.Cos(2 * Math.PI * fp * t));
            // Сигнал с частотой равной частоте заграждения
            var s_fs = new SamplesDigitalSignal(dt, samples_count, t => a0 * Math.Cos(2 * Math.PI * fs * t));
            // Сигнал с частотой равной половине частоты дискретизации
            var s_fd05 = new SamplesDigitalSignal(dt, samples_count, t => a0 * Math.Cos(2 * Math.PI * fd / 2 * t));

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
}