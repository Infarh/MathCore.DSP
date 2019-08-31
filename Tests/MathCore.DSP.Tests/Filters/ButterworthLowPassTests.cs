using System;
using System.Linq;
using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Service;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class ButterworthLowPassTests : UnitTest
    {
        [TestMethod]
        public void CreationTest()
        {
            const double fp = 0.05;    // Гц // Граничная частота полосы пропускания
            const double fs = 0.15;    // Гц // Граничная частота полосы запирания

            Assert.IsTrue(fp < fs);

            const double fd = 0.5;    // Гц // Частота дискретизации
            const double dt = 1 / fd; // 2с // Период дискретизации

            Assert.IsTrue(fp < fd / 2);

            //const double wp = Consts.pi2 * fp * dt; // 0.628318530717959 рад/с
            //const double ws = Consts.pi2 * fs * dt; // 1.884955592153876 рад/с

            const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 30; // Неравномерность в полосе пропускания (дБ)

            var eps_p = Math.Sqrt(Math.Pow(10, Rp / 10) - 1);
            var eps_s = Math.Sqrt(Math.Pow(10, Rs / 10) - 1);

            Assert.That.Value(eps_p).AreEqual(0.508847139909588, 4.441e-16);
            Assert.That.Value(eps_s).AreEqual(31.606961258558215);

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            Assert.That.Value(Gp).AreEqual(0.891250938133746, 4.45e-16);
            Assert.That.Value(Gs).AreEqual(0.0316227766016838, 4.45e-16);

            Assert.That.Value(0.891250938133746.In_dB()).AreEqual(-1, 4.67e-15);
            Assert.That.Value(0.0316227766016838.In_dB()).AreEqual(-30);

            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);  // Частота подавления аналогового протипа

            Assert.That.Value(Fp).AreEqual(0.051712575763384, 2.5e-16);
            Assert.That.Value(Fs).AreEqual(0.219057986225303, 3.89e-16);

            var Wp = Consts.pi2 * Fp;
            var Ws = Consts.pi2 * Fs;

            var k_eps = eps_s / eps_p;
            var k_W = Ws / Wp;

            Assert.That.Value(k_eps).GreaterThen(0);
            Assert.That.Value(k_W).GreaterThen(0);

            var double_n = Math.Log(k_eps) / Math.Log(k_W);
            // Порядок фильтра
            var N = (int)double_n;
            if (double_n - N > 0) N += 1;

            Assert.That.Value(N).AreEqual(3);

            var L = N / 2;
            var r = N % 2;

            Assert.That.Value(L).AreEqual(1);
            Assert.That.Value(r).AreEqual(1);
            Assert.That.Value(2 * L + r).AreEqual(N);

            var alpha = Math.Pow(eps_p, -1d / N);
            Assert.That.Value(alpha).AreEqual(1.252576388181026, 4.45e-16);

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

            Assert.That.Value(poles[0].Re).AreEqual(-1.252576388181026, 4.45e-16);
            Assert.That.Value(poles[0].Im).AreEqual(0);

            Assert.That.Value(poles[1].Re).AreEqual(-0.626288194090513, 4.45e-16);
            Assert.That.Value(poles[1].Im).AreEqual(1.084762972345327, 4.45e-16);

            Assert.That.Value(poles[2].Re).AreEqual(-0.626288194090513, 4.45e-16);
            Assert.That.Value(poles[2].Im).AreEqual(-1.084762972345327, 4.45e-16);

            var translated_poles = poles.Select(p => p * Wp).ToArray();

            Assert.That.Value(translated_poles[0].Re).AreEqual(-0.40698673955629);
            Assert.That.Value(translated_poles[0].Im).AreEqual(0);

            Assert.That.Value(translated_poles[1].Re).AreEqual(-0.203493369778145, 1.12e-16);
            Assert.That.Value(translated_poles[1].Im).AreEqual(0.352460855459148, 4.45e-16);

            Assert.That.Value(translated_poles[2].Re).AreEqual(-0.203493369778145, 1.12e-16);
            Assert.That.Value(translated_poles[2].Im).AreEqual(-0.352460855459148, 4.45e-16);

            var z_poles = translated_poles.Select(p => DigitalFilter.ToZ(p, dt)).ToArray();

            Assert.That.Value(z_poles[0].Re).AreEqual(0.421477504919999, 2.23e-16);
            Assert.That.Value(z_poles[0].Im).AreEqual(0);

            Assert.That.Value(z_poles[1].Re).AreEqual(0.530553579281761, 2.23e-16);
            Assert.That.Value(z_poles[1].Im).AreEqual(0.4482452811345, 4.45e-16);

            Assert.That.Value(z_poles[2].Re).AreEqual(0.530553579281761, 2.23e-16);
            Assert.That.Value(z_poles[2].Im).AreEqual(-0.4482452811345, 4.45e-16);

            Assert.That.Value(z_poles[0].Abs).LessThen(1);
            Assert.That.Value(z_poles[1].Abs).LessThen(1);
            Assert.That.Value(z_poles[2].Abs).LessThen(1);

            var kz = DigitalFilter.GetNomalizeCoefficient(translated_poles, dt);
            Assert.That.Value(kz).AreEqual(0.451944218734017, 1.12e-16);

            var WpN = Math.Pow(Wp, N);
            Assert.That.Value(WpN).AreEqual(0.034302685030762, 2.78e-16);

            var k = WpN * kz / eps_p;
            Assert.That.Value(k).AreEqual(0.030466713814018, 4.1e-16);

            var b = new double[N + 1];
            for (var i = 0; i < b.Length; i++)
                b[i] = k * SpecialFunctions.BinomialCoefficient(N, i);

            Assert.That.Value(b[0]).AreEqual(1 * k);
            Assert.That.Value(b[1]).AreEqual(3 * k);
            Assert.That.Value(b[2]).AreEqual(3 * k);
            Assert.That.Value(b[3]).AreEqual(1 * k);

            var a_complex = Polynom.Array.GetCoefficientsInverted(z_poles);
            Assert.That.Value(a_complex.Length).AreEqual(4);

            Assert.That.Value(a_complex[0].Im).AreEqual(0);
            Assert.That.Value(a_complex[1].Im).AreEqual(0);
            Assert.That.Value(a_complex[2].Im).AreEqual(0, 2.78e-1);
            Assert.That.Value(a_complex[3].Im).AreEqual(0);

            var a = a_complex.ToRe();

            Assert.That.Value(a[0]).AreEqual(1);
            Assert.That.Value(a[1]).AreEqual(-1.482584663483521, 6.67e-16);
            Assert.That.Value(a[2]).AreEqual(0.929643730192138, 1.12e-16);
            Assert.That.Value(a[3]).AreEqual(-0.203325356196476, 3.1e-16);

            var filter = new ButterworthLowPass(fp, fs, dt, Gp, Gs);

            var A = filter.A;
            var B = filter.B;

            Assert.That.Value(A.Count).AreEqual(a.Length);
            Assert.That.Value(B.Count).AreEqual(b.Length);

            Assert.That.Value(A[0]).AreEqual(a[0], 5.56e-17);
            Assert.That.Value(A[1]).AreEqual(a[1], 2.23e-16);
            Assert.That.Value(A[2]).AreEqual(a[2], 5.56e-17);
            Assert.That.Value(A[3]).AreEqual(a[3], 5.56e-17);

            Assert.That.Value(B[0]).AreEqual(b[0], 1.39e-17);
            Assert.That.Value(B[1]).AreEqual(b[1], 2.78e-17);
            Assert.That.Value(B[2]).AreEqual(b[2], 2.78e-17);
            Assert.That.Value(B[3]).AreEqual(b[3], 1.39e-17);
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

            Assert.That.Value(transmission_0.Abs).AreEqual(1, 2.23e-16);
            Assert.That.Value(transmission_fp.Abs).GreaterOrEqualsThen(Gp);
            Assert.That.Value(transmission_fs.Abs).LessOrEqualsThen(Gs);
            Assert.That.Value(transmission_fd05.Abs).AreEqual(0);
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
                0.030466713814017,
                0.136569624089457,
                0.265552982138031,
                0.303405053082752,
                0.230722080773791,
                0.11400006758795,
                0.016215356906167,
                -0.035027159443255,
                -0.043826129927865,
                -0.029116175709228,
                -0.009546418330117,
                0.004003333309531,
                0.008889991716749,
                0.007517462759302,
                0.003694729105447,
                0.000296747318716
            };
            var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();
            Assert.That.Value(impulse_response.Length).AreEqual(expected_impulse_response.Length);
            CollectionAssert.AreEqual(expected_impulse_response, impulse_response, GetComparer(8.44e-10));
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

            const int count = 1024;
            var s_0 = new SamplesDigitalSignal(dt, Enumerable.Repeat(1d, count));
            const double a0 = Consts.sqrt_2;
            var s_fp = new SamplesDigitalSignal(dt, count, t => a0 * Math.Cos(2 * Math.PI * fp * t));
            var s_fs = new SamplesDigitalSignal(dt, count, t => a0 * Math.Cos(2 * Math.PI * fs * t));
            var s_fd05 = new SamplesDigitalSignal(dt, count, t => a0 * Math.Cos(2 * Math.PI * fd / 2 * t));

            var y_0 = filter.ProcessIndividual(s_0);
            var y_fp = filter.ProcessIndividual(s_fp);
            var y_fs = filter.ProcessIndividual(s_fs);
            var y_fd05 = filter.ProcessIndividual(s_fd05);

            Assert.That.Value(y_0.Power).AreEqual(s_0.Power, 2.81e-3);
            Assert.That.Value(y_fp.Power).AreEqual(s_0.Power * Gp * Gp, 3.2e-3);
            Assert.That.Value(y_fs.Power).AreEqual(s_fs.Power * Gs * Gs, 1.84e-4);
            Assert.That.Value(y_fd05.Power).AreEqual(0, 1.33e-4);
        }
    }
}