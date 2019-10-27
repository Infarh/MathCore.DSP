using System;
using System.Linq;
using MathCore.DSP.Filters;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class EllipticLowPassTests : UnitTest
    {
        [TestMethod]
        public void Ctreation_Test()
        {
            const double fp = 100 / Math.PI / 2;
            const double fs = 150 / Math.PI / 2;
            const double fd = 100;
            const double dt = 1 / fd;

            const double wp = 2 * Math.PI * fp / fd; // 1
            const double ws = 2 * Math.PI * fs / fd; // 1.5

            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);
            Assert.That.Value(Fp).IsEqual(17.389348336409842, 3.56e-15);
            Assert.That.Value(Fs).IsEqual(29.653636313402, 3.56e-15);

            const double Rp = 1;
            const double Rs = 45;

            var eps_p = Math.Sqrt(Math.Pow(10, Rp / 10) - 1);
            var eps_s = Math.Sqrt(Math.Pow(10, Rs / 10) - 1);
            Assert.That.Value(eps_p).IsEqual(0.508847139909588, 4.45e-16);
            Assert.That.Value(eps_s).IsEqual(177.82512927503748, 4.45e-16);

            var k_eps = eps_s / eps_p;
            var k_W = Fs / Fp;
            Assert.That.Value(k_eps).IsEqual(349.46669702542425);
            Assert.That.Value(k_W).IsEqual(1.705275881518411, 2.23e-16);

            var k1 = eps_p / eps_s;
            var k = fp / fs;
            Assert.That.Value(k1).IsEqual(0.002861502994453, 4.28e-16);
            Assert.That.Value(k).IsEqual(0.666666666666667, 3.34e-16);

            var Kk = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(k);
            var Ksk = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(k);
            var Kk1 = SpecialFunctions.EllipticJacobi.FullEllipticIntegral(k1);
            var Ksk1 = SpecialFunctions.EllipticJacobi.FullEllipticIntegralComplimentary(k1);

            Assert.That.Value(Kk).IsEqual(1.809667495486589, 4.45e-16);
            Assert.That.Value(Ksk).IsEqual(1.904241416944999, 2.23e-16);
            Assert.That.Value(Kk1).IsEqual(1.570799542308087, 2.23e-16);
            Assert.That.Value(Ksk1).IsEqual(7.242715409944309, 8.89e-16);

            var double_N = Ksk1 * Kk / Kk1 / Ksk;
            Assert.That.Value(double_N).IsEqual(4.381849263936846);

            var N = (int)double_N;
            if (double_N - N > 0) N++;
            Assert.That.Value(N).IsEqual(5);

            var L = N / 2;
            var r = N % 2;
            Assert.That.Value(L).IsEqual(2);
            Assert.That.Value(r).IsEqual(1);

            double U(int i) => (2 * i - 1d) / N;
            var u = new double[L];
            for (var i = 0; i < L; i++)
                u[i] = U(i + 1);
            Assert.That.Value(u.Length).IsEqual(2);
            Assert.That.Value(u[0]).IsEqual(0.2);
            Assert.That.Value(u[1]).IsEqual(0.6);

            var k1p = Math.Sqrt(1 - k1 * k1);
            Assert.That.Value(k1p).IsEqual(0.999995905891925, 4.46e-16);

            var kp = k1p.Power(N) * u.Aggregate(1d, (P, ui) => P * SpecialFunctions.EllipticJacobi.sn_uk(ui, k1p).Power(4));
            Assert.That.Value(kp).IsEqual(0.641933634502708, 1.12e-16);

            k = Math.Sqrt(1 - kp * kp);
            Assert.That.Value(k).IsEqual(0.766760202993181, 4.45e-16);

            var im_pz = new double[L];
            for (var i = 0; i < L; i++)
                im_pz[i] = 1 / (k * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k));
            Assert.That.Value(im_pz.Length).IsEqual(2);
            Assert.That.Value(im_pz[0]).IsEqual(1.346819766817674, 4.45e-16);
            Assert.That.Value(im_pz[1]).IsEqual(1.945521905603307, 2.23e-16);

            var v0_complex = SpecialFunctions.EllipticJacobi.sn_inverse(new Complex(0, 1 / eps_p), k1) / N;
            Assert.That.Value(v0_complex.Im).IsEqual(0.181814340149935, 1.12e-16);
            Assert.That.Value(v0_complex.Re).IsEqual(0);

            var Pp = new Complex[N];
            var P0 = new Complex[N - r];

            if (r != 0) Pp[0] = Complex.i * SpecialFunctions.EllipticJacobi.sn_uk(v0_complex, k);
            for (var i = 0; i < L; i++)
            {
                var (p_im, p_re) = SpecialFunctions.EllipticJacobi.cd_uk(u[i] - v0_complex, k);

                Pp[r + 2 * i] = new Complex(-p_re, p_im);
                Pp[r + 2 * i + 1] = new Complex(-p_re, -p_im);

                var p0_im = 1 / (k * SpecialFunctions.EllipticJacobi.cd_uk(u[i], k));
                P0[2 * i] = new Complex(0, p0_im);
                P0[2 * i + 1] = new Complex(0, -p0_im);
            }

            Assert.That.Value(Pp.Length).IsEqual(5);
            Assert.That.Value(P0.Length).IsEqual(4);

            Assert.That.Value(Pp[0].Re).IsEqual(-0.364129349944568, 5.56e-17);
            Assert.That.Value(Pp[0].Im).IsEqual(0);

            Assert.That.Value(Pp[1].Re).IsEqual(-0.056735988486377, 3.62e-16);
            Assert.That.Value(Pp[1].Im).IsEqual(0.997076970459406, 5.56e-16);
            Assert.That.Value(Pp[2].Re).IsEqual(-0.056735988486377, 3.62e-16);
            Assert.That.Value(Pp[2].Im).IsEqual(-0.997076970459406, 5.56e-16);

            Assert.That.Value(Pp[3].Re).IsEqual(-0.223929681527417, 3.34e-16);
            Assert.That.Value(Pp[3].Im).IsEqual(0.715628507541036, 2.23e-16);
            Assert.That.Value(Pp[4].Re).IsEqual(-0.223929681527417, 3.34e-16);
            Assert.That.Value(Pp[4].Im).IsEqual(-0.715628507541036, 3.34e-16);

            Assert.That.Value(P0[0].Re).IsEqual(0);
            Assert.That.Value(P0[0].Im).IsEqual(1.346819766817674, 4.45e-16);

            Assert.That.Value(P0[1].Re).IsEqual(0);
            Assert.That.Value(P0[1].Im).IsEqual(-1.346819766817674, 4.45e-16);

            var numirator_coefficients = Polynom.Array.GetCoefficients(P0);
            var denomirator_coefficients = Polynom.Array.GetCoefficients(Pp);

            Assert.That.Value(numirator_coefficients.Length).IsEqual(5);
            Assert.That.Value(denomirator_coefficients.Length).IsEqual(6);

            var (B, numirator_coefficients_im) = numirator_coefficients;
            var (A, denomirator_coefficients_im) = denomirator_coefficients;

            CollectionAssert.That.Collection(numirator_coefficients_im).ElementsAreEqualTo(0);
            CollectionAssert.That.Collection(denomirator_coefficients_im).ElementsAreEqualTo(0, 5.56e-17);

            CollectionAssert.AreEqual(new[]
            {
                6.865801033915977,
                0,
                5.598978969473136,
                0,
                1
            }, B, GetComparer(5.33e-15));

            CollectionAssert.AreEqual(new[]
            {
                0.204202406237771,
                0.746680133688232,
                1.096907612426762,
                1.814866823763763,
                0.925460689972156,
                1
            }, A, GetComparer(6.67e-16));

            var norm_k = B[0] / A[0];
            Assert.That.Value(norm_k).IsEqual(33.62252757159741, 2.85e-14);

            B = B.Select(v => v / norm_k).ToArray();
            CollectionAssert.That.Collection(B).IsEqualTo(new[]
            {
                0.204202406237771,
                0,
                0.166524630176908,
                0,
                0.029741963862489
            }, 5.1e-16);
        }
    }
}
