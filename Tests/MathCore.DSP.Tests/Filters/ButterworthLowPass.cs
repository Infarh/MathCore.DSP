using System;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Service;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class ButterworthLowPass : UnitTest
    {
        [TestMethod]
        public void Creation()
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

            Assert.That.Value(eps_p).IsEqual(0.50884713990958752);
            Assert.That.Value(eps_s).IsEqual(31.606961258558215);

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            Assert.That.Value(Gp).IsEqual(0.89125093813374556);
            Assert.That.Value(Gs).IsEqual(0.031622776601683791);

            Assert.That.Value(Gp.In_dB()).IsEqual(-1, 4.67e-15);
            Assert.That.Value(Gs.In_dB()).IsEqual(-30);

            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);  // Частота подавления аналогового прототипа

            Assert.That.Value(Fp).IsEqual(0.051712575763384123);
            Assert.That.Value(Fs).IsEqual(0.2190579862253032);

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

            var alpha = eps_p.Pow(-1d / N);
            Assert.That.Value(alpha).IsEqual(1.2525763881810263);

            var th0 = Math.PI / N;

            var poles = new Complex[N];
            if (r != 0) poles[0] = -alpha;
            for (var i = r; i < poles.Length; i += 2)
            {
                var w = th0 * (i + 1 - r - 0.5);
                poles[i] = (-alpha * Math.Sin(w), alpha * Math.Cos(w));
                poles[i + 1] = poles[i].ComplexConjugate;
            }

            Assert.That.Collection(poles).IsEqualTo(
                -1.2525763881810263,
                (-0.62628819409051306, 1.0847629723453271),
                (-0.62628819409051306, -1.0847629723453271)
                );

            var translated_poles = poles.ToArray(p => p * Wp)
               .AssertThatCollection()
               .IsEqualTo(
                    -0.40698673955629,
                    (-0.20349336977814494, 0.35246085545914824),
                    (-0.20349336977814494, -0.35246085545914824))
               .ActualValue;

            var z_poles = translated_poles.ToArray(p => DigitalFilter.ToZ(p, dt))
               .AssertThatCollection()
               .IsEqualTo(
                     0.42147750491999925,
                     (0.53055357928176117, 0.44824528113449957),
                     (0.53055357928176117, -0.44824528113449957))
               .AllItems(z_pole => z_pole.Where(z => z.Abs < 1))
               .ActualValue;

            var kz = DigitalFilter.GetNormalizeCoefficient(translated_poles, dt)
                .AssertThanValue()
                .IsEqual(0.45194421873401691)
                .ActualValue;

            var WpN = Wp.Pow(N)
               .AssertThanValue()
               .IsEqual(0.034302685030761962)
               .ActualValue;

            var k = WpN * kz / eps_p;
            Assert.That.Value(k).IsEqual(0.030466713814017589);

            var b = Enumerable.Range(0, N + 1).ToArray(i => k * SpecialFunctions.BinomialCoefficient(N, i));

            Assert.That.Collection(b).IsEqualTo(new[]
            {
                1 * k,
                3 * k,
                3 * k,
                1 * k
            });

            var a_complex = Polynom.Array.GetCoefficientsInverted((Complex[])z_poles);
            Assert.That.Collection(a_complex)
               .CountEquals(4)
               .AllItems(a_value => a_value.Where(z => z.Im).IsEqual(0, 1e-16));

            var a = a_complex.ToRe() ?? throw new AssertFailedException("Отсутствует ссылка на массив вещественных значений");

            Assert.That.Collection(a).ValuesAreEqual(
                1,
                -1.4825846634835216,
                0.92964373019213786,
                -0.20332535619647568
            );

            var filter = new DSP.Filters.ButterworthLowPass(fp, fs, dt, Gp, Gs);

            Assert.That.Collection(filter.A).IsEqualTo(a, 1e-15);
            Assert.That.Collection(filter.B).IsEqualTo(b, 1e-16);
        }

        [TestMethod]
        public void TransmissionCoefficient()
        {
            const double fp = 0.1;    // Гц // Граничная частота полосы пропускания
            const double fs = 0.3;    // Гц // Граничная частота полосы запирания
            const double fd = 1;      // Гц // Частота дискретизации
            const double dt = 1 / fd; // 2с // Период дискретизации
            const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 30; // Неравномерность в полосе пропускания (дБ)
            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new DSP.Filters.ButterworthLowPass(fp, fs, dt, Gp, Gs);

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            Assert.That.Value(transmission_0.Abs).IsEqual(1, 2.23e-16);
            Assert.That.Value(transmission_fp.Abs).IsEqual(Gp, 6.67e-016);
            Assert.That.Value(transmission_fs.Abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05.Abs).IsEqual(0, 4.27e-34);
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

            var filter = new DSP.Filters.ButterworthLowPass(fp, fs, dt, Gp, Gs);

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
        public void SignallProcessing()
        {
            const double fp = 0.1;    // Гц // Граничная частота полосы пропускания
            const double fs = 0.3;    // Гц // Граничная частота полосы запирания
            const double fd = 1;      // Гц // Частота дискретизации
            const double dt = 1 / fd; // 2с // Период дискретизации
            const double Rp = 1;      // Неравномерность в полосе пропускания (дБ)
            const double Rs = 30;     // Неравномерность в полосе пропускания (дБ)
            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new DSP.Filters.ButterworthLowPass(fp, fs, dt, Gp, Gs);

            const int samples_count = 1024;
            // Сигнал s0(t) = 1
            var s_0 = new SamplesDigitalSignal(dt, Enumerable.Repeat(1d, samples_count));

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
}