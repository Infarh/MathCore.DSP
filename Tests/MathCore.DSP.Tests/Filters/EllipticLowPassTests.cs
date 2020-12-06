using System;
using System.Collections.ObjectModel;
using System.Linq;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class EllipticLowPassTests : UnitTest
    {
        [TestMethod]
        public void Creation_Test()
        {
            //https://ru.dsplib.org/content/filter_ellip_ap/filter_ellip_ap.html

            const double pi2 = 2 * Math.PI;

            const double fd = 5000;
            const double dt = 1 / fd;

            const double fp = fd / pi2;
            const double fs = 1.5 * fp;

            //const double wp = 2 * Math.PI * fp / fd; // 1
            //const double ws = 2 * Math.PI * fs / fd; // 1.5

            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);
            Assert.That.Value(Fp).IsEqual(869.4674168204921, 3.56e-15);
            Assert.That.Value(Fs).IsEqual(1482.6818156701001, 3.56e-15);

            var Wp = Consts.pi2 * Fp;
            var Ws = Consts.pi2 * Fs;

            const double Rp = 1;
            const double Rs = 45;

            var eps_p = Math.Sqrt(Math.Pow(10, Rp / 10) - 1);
            var eps_s = Math.Sqrt(Math.Pow(10, Rs / 10) - 1);
            Assert.That.Value(eps_p).IsEqual(0.508847139909588, 4.45e-16);
            Assert.That.Value(eps_s).IsEqual(177.82512927503748, 4.45e-16);

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            //var k_eps = eps_s / eps_p;
            //var k_W = Fs / Fp;
            //Assert.That.Value(k_eps).IsEqual(349.46669702542425);
            //Assert.That.Value(k_W).IsEqual(1.705275881518411, 2.23e-16);

            var k_W = fp / fs;
            var k_eps = eps_p / eps_s;
            Assert.That.Value(k_W).IsEqual(0.666666666666667, 3.34e-16);
            Assert.That.Value(k_eps).IsEqual(0.002861502994453, 4.28e-16);

            var K_w = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(k_W);
            var T_w = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(k_W);
            var K_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(k_eps);
            var T_eps = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(k_eps);

            Assert.That.Value(K_w).IsEqual(1.809667495486589, 4.45e-16);
            Assert.That.Value(T_w).IsEqual(1.904241416944999, 2.23e-16);
            Assert.That.Value(K_eps).IsEqual(1.570799542308087, 2.23e-16);
            Assert.That.Value(T_eps).IsEqual(7.242715409944309, 8.89e-16);

            // Оценка снизу порядка фильтра
            var double_N = T_eps * K_w / K_eps / T_w;
            Assert.That.Value(double_N).IsEqual(4.381849263936846);

            var N = (int)Math.Ceiling(double_N); // Порядок фильтра
            Assert.That.Value(N).IsEqual(5);

            var L = N / 2;
            var r = N % 2;
            Assert.That.Value(L).IsEqual(2); // Число чётных (парных) полюсов
            Assert.That.Value(r).IsEqual(1); // Есть один нечётный полюс

            // Эллиптический модуль
            double U(int i) => (2 * i - 1d) / N;
            var u = new double[L];
            for (var i = 0; i < L; i++)
                u[i] = U(i + 1);

            Assert.That.Collection(u).ValuesAreEqual(0.2, 0.6);

            var m = (1 - k_eps * k_eps).Sqrt();
            Assert.That.Value(m).IsEqual(0.999995905891925, 4.46e-16);

            
            var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * SpecialFunctions.EllipticJacobi.sn_uk(ui, m).Power(4));
            Assert.That.Value(kp).IsEqual(0.641933634502708, 1.12e-16);

            k_W = (1 - kp * kp).Sqrt();
            Assert.That.Value(k_W).IsEqual(0.766760202993181, 4.45e-16);

            var im_pz = new double[L];
            for (var i = 0; i < L; i++)
                im_pz[i] = 1 / (k_W * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k_W));
            Assert.That.Collection(im_pz)
               .ValuesAreEqual(
                    1.3468197668176745,
                    1.9455219056033073);
            
            var v0_complex = SpecialFunctions.EllipticJacobi.sn_inverse(new Complex(0, 1 / eps_p), k_eps) / N;
            Assert.That.Value(v0_complex.Im).IsEqual(0.181814340149935, 1.12e-16);
            Assert.That.Value(v0_complex.Re).IsEqual(0);

            var poles = new Complex[N];
            var zeros = new Complex[N - r];

            if (r != 0) poles[0] = Complex.i * SpecialFunctions.EllipticJacobi.sn_uk(v0_complex, k_W);
            for (var i = 0; i < L; i++)
            {
                var (p_im, p_re) = SpecialFunctions.EllipticJacobi.cd_uk(u[i] - v0_complex, k_W);

                poles[r + 2 * i] = new Complex(-p_re, p_im);
                poles[r + 2 * i + 1] = new Complex(-p_re, -p_im);

                var p0_im = 1 / (k_W * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k_W));
                zeros[2 * i] = new Complex(0, p0_im);
                zeros[2 * i + 1] = new Complex(0, -p0_im);
            }

            Assert.That.Value(poles.Length).IsEqual(5);
            Assert.That.Value(zeros.Length).IsEqual(4);

            // Полюса
            Assert.That.Collection(poles).IsEqualTo(
                new Complex(-0.36412934994456797),
                new Complex(-0.05673598848637736, +0.9970769704594055),
                new Complex(-0.05673598848637736, -0.9970769704594055),
                new Complex(-0.2239296815274167, +0.7156285075410361),
                new Complex(-0.2239296815274167, -0.7156285075410361)
                );

            // Нули
            Assert.That.Collection(zeros).IsEqualTo(
                new Complex(0, +1.3468197668176745),
                new Complex(0, -1.3468197668176745),
                new Complex(0, +1.9455219056033073),
                new Complex(0, -1.9455219056033073)
                );

            var analog_numerator_coefficients = Polynom.Array.GetCoefficientsInverted(zeros);
            var analog_denominator_coefficients = Polynom.Array.GetCoefficientsInverted(poles);

            Assert.That.Value(analog_numerator_coefficients.Length).IsEqual(5);
            Assert.That.Value(analog_denominator_coefficients.Length).IsEqual(6);

            var (analog_b, numerator_coefficients_im) = analog_numerator_coefficients;
            var (analog_a, denominator_coefficients_im) = analog_denominator_coefficients;

            CollectionAssert.That.Collection(numerator_coefficients_im).ElementsAreEqualTo(0);
            CollectionAssert.That.Collection(denominator_coefficients_im).ElementsAreEqualTo(0, 5.56e-17);

            Assert.That.Collection(analog_b).ValuesAreEqual(
                1,
                0,
                5.598978969473139,
                0,
                6.865801033915982);

            Assert.That.Collection(analog_a).ValuesAreEqual(
                1,
                0.925460689972156,
                1.8148668237637633,
                1.0969076124267614,
                0.7466801336882324,
                0.2042024062377706);

            var norm_k = analog_b[^1] / analog_a[^1];
            Assert.That.Value(norm_k).IsEqual(33.62252757159741, 2.85e-14);

            analog_b = analog_b.ToArray(v => v / norm_k);
            Assert.That.Collection(analog_b).ValuesAreEqual(
                0.029741963862489267,
                0,
                0.1665246301769075,
                0,
                0.20420240623777058);

            var translated_zeros = poles.ToArray(p => p * Wp);
            var translated_poles = poles.ToArray(p => p * Wp);

            var kz = DigitalFilter.GetNormalizeCoefficient(translated_poles, dt);

            var WpN = Wp.Pow(N);

            var k = WpN * kz / eps_p;

            var z_zeros = translated_poles.ToArray(p => DigitalFilter.ToZ(p, dt));
            var z_poles = translated_poles.ToArray(p => DigitalFilter.ToZ(p, dt));

            var b_complex = Polynom.Array.GetCoefficientsInverted(z_poles);
            var a_complex = Polynom.Array.GetCoefficientsInverted(z_poles);

            var b = b_complex.ToRe();
            var a = a_complex.ToRe();

            var H0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(a, b, 0);
            var Hp = DoubleArrayDSPExtensions.GetTransmissionCoefficient(a, b, fp);
            var Hs = DoubleArrayDSPExtensions.GetTransmissionCoefficient(a, b, fs);
            var Hd2 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(a, b, fd / 2);
        }

        [TestMethod, Ignore]
        public void TransmissionCoefficientTest()
        {
            const double pi2 = 2 * Math.PI;
            //const double fd = 5000;
            //const double dt = 1 / fd;

            //const double fp = fd / pi2;
            //const double fs = 1.5 * fp;

            const double fd = 5000;      // Гц // Частота дискретизации
            const double dt = 1 / fd;    // с  // Период дискретизации
            const double fp = fd / pi2;  // Гц // Граничная частота полосы пропускания
            const double fs = 1.5 * fp;  // Гц // Граничная частота полосы запирания
            const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 30; // Неравномерность в полосе пропускания (дБ)

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new EllipticLowPassFilter(fp, fs, dt, Gp, Gs);

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            Assert.That.Value(transmission_0.Abs).IsEqual(1, 2.23e-16);
            Assert.That.Value(transmission_fp.Abs).GreaterOrEqualsThan(Gp);
            Assert.That.Value(transmission_fs.Abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05.Abs).IsEqual(0, 4.27e-34);
        }
    }
}
