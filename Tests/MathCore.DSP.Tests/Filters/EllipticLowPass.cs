using System;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Service;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions.EllipticJacobi;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class EllipticLowPass : UnitTest
    {
        [TestMethod]
        public void Creation()
        {
            //https://ru.dsplib.org/content/filter_ellip_ap/filter_ellip_ap.html

            const double pi2 = 2 * Math.PI;

            const double fd = 5000;         // Частота дискретизации
            const double dt = 1 / fd;       // Период дискретизации

            const double fp = fd / pi2;     // Граничная частота конца интервала пропускания
            const double fs = 1.5 * fp;     // Граничная частота начала интервала подавления

            const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
            const double Rs = 45;   // Уровень подавления более 45 дБ

            var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
            var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

            #region Аналитический расчёт фильтра

            //const double wp = 2 * Math.PI * fp / fd; // 1
            //const double ws = 2 * Math.PI * fs / fd; // 1.5

            // Рассчитываем частоты цифрового фильтра
            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt).AssertThanValue().IsEqual(869.46741682049208).ActualValue;
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt).AssertThanValue().IsEqual(1482.6818156701001).ActualValue;

            // Круговые частоты
            var Wp = Consts.pi2 * Fp;
            //var Ws = Consts.pi2 * Fs;

            // Допуск на АЧХ в интервале пропускания
            var eps_p = (Math.Pow(10, Rp / 10) - 1).Sqrt().AssertThanValue().IsEqual(0.50884713990958752).ActualValue;
            // Допуск на АЧХ в интервале подавления
            var eps_s = (Math.Pow(10, Rs / 10) - 1).Sqrt().AssertThanValue().IsEqual(177.82512927503748).ActualValue;

            //var k_eps = eps_s / eps_p;
            //var k_W = Fs / Fp;
            //Assert.That.Value(k_eps).IsEqual(349.46669702542425);
            //Assert.That.Value(k_W).IsEqual(1.705275881518411, 2.23e-16);

            var k_W = fp / fs;
            var k_eps = eps_p / eps_s;
            Assert.That.Value(k_W).IsEqual(0.66666666666666663);
            Assert.That.Value(k_eps).IsEqual(0.0028615029944534269);

            var K_w = FullEllipticIntegral(k_W)
               .AssertThanValue().IsEqual(1.8096674954865886)
               .ActualValue;

            var T_w = FullEllipticIntegralComplimentary(k_W)
               .AssertThanValue().IsEqual(1.9042414169449993)
               .ActualValue;

            var K_eps = FullEllipticIntegral(k_eps)
               .AssertThanValue().IsEqual(1.5707995423080867)
               .ActualValue;

            var T_eps = FullEllipticIntegralComplimentary(k_eps)
               .AssertThanValue().IsEqual(7.2427154099443083)
               .ActualValue;

            // Оценка снизу порядка фильтра
            var double_N = T_eps * K_w / (K_eps * T_w);
            Assert.That.Value(double_N).IsEqual(4.381849263936846);

            var N = (int)Math.Ceiling(double_N); // Порядок фильтра
            Assert.That.Value(N).IsEqual(5);

            var L = N / 2;  // Число комплексно сопряжённых полюсов
            var r = N % 2;  // Число (0 или 1) действительных полюсов - (чётность фильтра)
            Assert.That.Value(L).IsEqual(2); // Число чётных (парных) полюсов
            Assert.That.Value(r).IsEqual(1); // Есть один нечётный полюс

            // Эллиптический модуль
            var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);

            Assert.That.Collection(u).ValuesAreEqual(0.2, 0.6);

            var m = (1 - k_eps * k_eps).Sqrt().AssertThanValue().IsEqual(0.99999590589192544).ActualValue;

            var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));
            Assert.That.Value(kp).IsEqual(0.64193363450270813);

            k_W = (1 - kp * kp).Sqrt().AssertThanValue().IsEqual(0.7667602029931806).ActualValue;

            var im_pz = Enumerable.Range(0, L).ToArray(i => 1 / (k_W * cd_uk(u[i], k_W)));

            Assert.That.Collection(im_pz)
               .ValuesAreEqual(
                    1.3468197668176745,
                    1.9455219056033073);

            var v0_complex = sn_inverse((0, 1 / eps_p), k_eps) / N;
            Assert.That.Value(v0_complex).IsEqual((0, 0.18181434014993489));

            var zeros = new Complex[N - r]; // Массив нулей (на r меньше числа полюсов)
            var poles = new Complex[N];     // Массив полюсов

            // Если фильтр нечётный, то первым полюсом будет действительный полюс
            if (r != 0) poles[0] = Complex.i * sn_uk(v0_complex, k_W);
            for (var i = 0; i < L; i++)
            {
                // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
                var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_W);

                poles[r + 2 * i] = (-p_re, p_im);
                poles[r + 2 * i + 1] = poles[r + 2 * i].ComplexConjugate;

                var p0_im = 1 / (k_W * cd_uk(u[i], k_W));
                zeros[2 * i] = (0, p0_im);
                zeros[2 * i + 1] = zeros[2 * i].ComplexConjugate;
            }

            // Полюса
            Assert.That.Collection(poles).IsEqualTo(
                (-0.36412934994456797),
                (-0.05673598848637736, +0.9970769704594055),
                (-0.05673598848637736, -0.9970769704594055),
                (-0.2239296815274167, +0.7156285075410361),
                (-0.2239296815274167, -0.7156285075410361)
                );

            // Нули
            Assert.That.Collection(zeros).IsEqualTo(
                (0, +1.3468197668176745),
                (0, -1.3468197668176745),
                (0, +1.9455219056033073),
                (0, -1.9455219056033073)
                );


            // Рассчитываем коэффициенты полиномов числителя и знаменателя передаточной функции
            // аналогового прототипа H(p) = Q(p) / P(p)
            var analog_numerator_coefficients = GetCoefficientsInverted(zeros);     // Q(p)
            var analog_denominator_coefficients = GetCoefficientsInverted(poles);   // P(p)

            var (analog_b, numerator_coefficients_im) = analog_numerator_coefficients;
            var (analog_a, denominator_coefficients_im) = analog_denominator_coefficients;

            CollectionAssert.That.Collection(numerator_coefficients_im).ElementsAreEqualTo(0);
            CollectionAssert.That.Collection(denominator_coefficients_im).ElementsAreEqualTo(0, 5.56e-17);

            Assert.That.Collection(analog_b)
               .ValuesAreEqual(
                    1,
                    0,
                    5.598978969473139,
                    0,
                    6.865801033915982);

            Assert.That.Collection(analog_a)
               .ValuesAreEqual(
                    1,
                    0.925460689972156,
                    1.8148668237637633,
                    1.0969076124267614,
                    0.7466801336882324,
                    0.2042024062377706);

            var norm_k = analog_b[^1] / analog_a[^1];
            Assert.That.Value(norm_k).IsEqual(33.62252757159741, 2.85e-14);

            analog_b = analog_b.ToArray(v => v / norm_k);
            Assert.That.Collection(analog_b)
               .ValuesAreEqual(
                    0.029741963862489267,
                    0,
                    0.1665246301769075,
                    0,
                    0.20420240623777058);

            var translated_zeros = zeros.ToArray(p => p * Wp);
            var translated_poles = poles.ToArray(p => p * Wp);

            Assert.That.Collection(translated_zeros)
               .IsEqualTo(
                    (0, 7357.709919833289),
                    (0, -7357.709919833289),
                    (0, 10628.434610767226),
                    (0, -10628.434610767226)
                    );

            Assert.That.Collection(translated_poles)
               .IsEqualTo(
                    -1989.2477049991837,
                    (-309.95011773856584, 5447.0563152787672),
                    (-309.95011773856584, -5447.0563152787672),
                    (-1223.3334256835481, 3909.4963547286384),
                    (-1223.3334256835481, -3909.4963547286384));

            var z_zeros = translated_zeros.ToArray(z => DigitalFilter.ToZ(z, dt));
            var z_poles = translated_poles.ToArray(z => DigitalFilter.ToZ(z, dt));

            Assert.That.Collection(z_zeros)
               .IsEqualTo(
                    (0.29755628730678862, 0.95470427666592117),
                    (0.29755628730678862, -0.95470427666592117),
                    (-0.060872472663867229, 0.9981455515463598),
                    (-0.060872472663867229, -0.9981455515463598));

            Assert.That.Collection(z_poles)
               .IsEqualTo(
                    0.66816138027247152,
                    (0.516553916670814, 0.80124098515759445),
                    (0.516553916670814, -0.80124098515759445),
                    (0.58917408994119846, 0.553567293780718),
                    (0.58917408994119846, -0.553567293780718));

            if (r > 0)
            {
                Array.Resize(ref z_zeros, z_zeros.Length + 1);
                Array.Copy(z_zeros, 0, z_zeros, 1, z_zeros.Length - 1);
                z_zeros[0] = -1;
            }

            var G_norm = (r > 0 ? 1 : 1 / (1 + eps_p * eps_p).Sqrt())
                / (z_zeros.Aggregate(Complex.Real, (Z, z) => Z * (1 - z))
                    / z_poles.Aggregate(Complex.Real, (Z, z) => Z * (1 - z))).Abs;

            Assert.That.Value(G_norm).IsEqual(0.023163864337949005);

            var (B, im_b) = GetCoefficientsInverted(z_zeros).ToArray(b => b * G_norm);
            var (A, im_a) = GetCoefficientsInverted(z_poles);

            Assert.That.Collection(im_a).ElementsAreEqualTo(0, 1e-15);
            Assert.That.Collection(im_b).ElementsAreEqualTo(0, 1e-15);

            // Проверяем коэффициенты передачи рассчитанного фильтра
            var H0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0).Abs;
            var Hp = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fp / fd).Abs;
            var Hs = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fs / fd).Abs;
            var Hd = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0.5).Abs;

            const double eps = 1e-14;
            Assert.That.Value(H0).IsEqual(1, eps);                 // Коэффициент на нулевой частоте должен быть равен 1
            Assert.That.Value(Hp).GreaterOrEqualsThan(Gp, eps);    // Коэффициент передачи на граничной частоте конца интервала пропускания должен быть не меньше заданного значения Gp
            Assert.That.Value(Hs).LessOrEqualsThan(Gs);            // Коэффициент передачи на граничной частоте начала интервала подавления должен быть не больше заданного значения Gs
            Assert.That.Value(Hd).IsEqual(0, eps);                 // Коэффициент передачи на частоте, равной половине частоты дискретизации (соответствующей бесконечно большой частоте) должен быть равен нулю

            #endregion

            var filter = new DSP.Filters.EllipticLowPass(fp, fs, dt, Gp, Gs);

            Assert.That.Collection(filter.A).IsEqualTo(A, 1e-15);
            Assert.That.Collection(filter.B).IsEqualTo(B, 1e-16);
        }

        [TestMethod]
        public void TransmissionCoefficient()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;      // Гц // Частота дискретизации
            const double dt = 1 / fd;    // с  // Период дискретизации
            const double fp = fd / pi2;  // Гц // Граничная частота полосы пропускания
            const double fs = 1.5 * fp;  // Гц // Граничная частота полосы запирания
            const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 45; // Неравномерность в полосе пропускания (дБ)

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new DSP.Filters.EllipticLowPass(fp, fs, dt, Gp, Gs);

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            var transmission_0_abs = transmission_0.Abs;
            var transmission_fp_abs = transmission_fp.Abs;
            var transmission_fs_abs = transmission_fs.Abs;
            var transmission_fd05_abs = transmission_fd05.Abs;

            const double eps = 2.18e-14;
            Assert.That.Value(transmission_0_abs).IsEqual(1, eps);
            Assert.That.Value(transmission_fp_abs).IsEqual(Gp, eps);
            Assert.That.Value(transmission_fs_abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05_abs).IsEqual(0, eps);
        }

        [TestMethod]
        public void ImpulseResponse()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;      // Гц // Частота дискретизации
            const double dt = 1 / fd;    // с  // Период дискретизации
            const double fp = fd / pi2;  // Гц // Граничная частота полосы пропускания
            const double fs = 1.5 * fp;  // Гц // Граничная частота полосы запирания
            const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 45; // Неравномерность в полосе пропускания (дБ)
            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new DSP.Filters.EllipticLowPass(fp, fs, dt, Gp, Gs);

            double[] expected_impulse_response =
            {
                +0.02316402765395902,
                +0.07890235175222579,
                +0.16227581269333505,
                +0.2485315665700753 ,
                +0.2805442139435834 ,
                +0.22795965026314102,
                +0.10322811584198478,
                -0.03550289833716826,
                -0.11558449693089311,
                -0.09973469008670285,
                -0.01432462412676058,
                +0.07033262928146783,
                +0.09354381942841251,
                +0.0480481590295324 ,
                -0.02081136468031311,
                -0.05692503754813389
            };
            var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();

            const double eps = 9.4e-7;
            Assert.That.Collection(impulse_response)
               .IsEqualTo(expected_impulse_response, eps);
        }

        [TestMethod]
        public void SignallProcessing()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;      // Гц // Частота дискретизации
            const double dt = 1 / fd;    // с  // Период дискретизации
            const double fp = fd / pi2;  // Гц // Граничная частота полосы пропускания
            const double fs = 1.5 * fp;  // Гц // Граничная частота полосы запирания
            const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 45; // Неравномерность в полосе пропускания (дБ)
            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new DSP.Filters.EllipticLowPass(fp, fs, dt, Gp, Gs);

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
            Assert.That.Value(y_0.Power).IsEqual(s_0.Power, 3.94e-003);
            // На граничной частоте сигнал должен быть ослаблен на коэффициент Gp (на Gp^2 по мощности)
            Assert.That.Value(y_fp.Power).IsEqual(s_0.Power * Gp * Gp, 2.26e-2);
            // На частоте заграждения сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
            Assert.That.Value(y_fs.Power).LessOrEqualsThan(s_fs.Power * Gs * Gs, 2.34e-4);
            // На частоте в половину частоты дискретизации сигнал должен быть подавлен
            Assert.That.Value(y_fd05.Power).IsEqual(0, 1.57e-4);
        }
    }
}
