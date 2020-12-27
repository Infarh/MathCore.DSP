using System;
using System.Linq;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class ChecyshevLowPass : UnitTest
    {
        private static double arsh(double x) => Math.Log(x + Math.Sqrt(x * x + 1));
        private static double arch(double x) => Math.Log(x + Math.Sqrt(x * x - 1));

        [TestMethod]
        public void TypeI_Creation()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Частота дискретизации
            const double dt = 1 / fd;           // Период дискретизации

            const double fp = 1.0 * fd / pi2;   // Граничная частота конца интервала пропускания
            const double fs = 1.6 * fd / pi2;   // Граничная частота начала интервала подавления

            const double Rp = 1.5;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 35; // Неравномерность в полосе пропускания (дБ)

            #region Аналитический расчёт

            Assert.IsTrue(fp < fs);
            Assert.IsTrue(fp < fd / 2);

            //const double wp = Consts.pi2 * fp * dt; // 0.628318530717959 рад/с
            //const double ws = Consts.pi2 * fs * dt; // 1.884955592153876 рад/с

            var eps_p = (10d.Pow(Rp / 10) - 1).Sqrt();
            var eps_s = (10d.Pow(Rs / 10) - 1).Sqrt();

            Assert.That.Value(eps_p).IsEqual(0.6422908567173865);
            Assert.That.Value(eps_s).IsEqual(56.225240418946896);

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            Assert.That.Value(Gp).IsEqual(0.84139514164519513);
            Assert.That.Value(Gs).IsEqual(0.017782794100389229);

            Assert.That.Value(Gp.In_dB()).IsEqual(-1.5, 4.67e-15);
            Assert.That.Value(Gs.In_dB()).IsEqual(-35);

            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);  // Частота подавления аналогового прототипа

            Assert.That.Value(Fp).IsEqual(869.46741682049208);
            Assert.That.Value(Fs).IsEqual(1638.7206595257194);

            var Wp = Consts.pi2 * Fp;
            var Ws = Consts.pi2 * Fs;

            var k_eps = eps_s / eps_p;
            var k_W = Ws / Wp;

            Assert.That.Value(k_eps).IsEqual(87.5385969314623);
            Assert.That.Value(k_W).IsEqual(1.8847407364824178);

            var N = (int)Math.Ceiling(arch(k_eps) / arch(k_W)); // Порядок фильтра
            Assert.That.Value(N).IsEqual(5);

            var r = N % 2;                              // Нечётность порядка фильтра
            var dth = Math.PI / N;                      // Угловой шаг между полюсами
            var beta = arsh(1 / eps_p) / N;
            Assert.That.Value(beta).IsEqual(0.24518628509618212);

            var sh = Math.Sinh(beta);
            var ch = Math.Cosh(beta);
            var poles = new Complex[N];                 // Массив полюсов фильтра
            if (r != 0) poles[0] = -sh;                 // Если порядок фильтра нечётный, то первым добавляем центральный полюс
            for (var i = r; i < poles.Length; i += 2)   // Расчёт полюсов
            {
                var n = (i - r) / 2 + 1;
                var th = dth * (n - 0.5);

                var sin = Math.Sin(th);
                var cos = Math.Cos(th);
                poles[i] = new Complex(-sh * sin, ch * cos);
                poles[i + 1] = poles[i].ComplexConjugate;
            }

            Assert.That.Collection(poles).IsEqualTo(
                -0.24765029577598821,
                (-0.076528150056762612, +0.979787021976867),
                (-0.076528150056762612, -0.979787021976867),
                (-0.20035329794475673, +0.605541681317744),
                (-0.20035329794475673, -0.605541681317744)
                );

            var z_poles = DigitalFilter.ToZArray(poles, dt, Wp);

            Assert.That.Collection(z_poles).IsEqualTo(
                0.76166135868560936,
                (0.51881789309466986, +0.78033858156608826),
                (0.51881789309466986, -0.78033858156608826),
                (0.65550339923693535, +0.49362618841054584),
                (0.65550339923693535, -0.49362618841054584)
            );

            var (A, a_im) = Polynom.Array.GetCoefficientsInverted(z_poles);

            Assert.That.Collection(a_im).ElementsAreEqualTo(0, 2.3e-16);
            Assert.That.Collection(A).ValuesAreEqual(
                +1,
                -3.11030394334882,
                +4.700669700507134,
                -4.067694193782959,
                +2.000259228659479,
                -0.4503476466802788);

            var G_norm = (r > 0 ? 1 : Gp)
                / (2.Power(N) / z_poles.Aggregate(Complex.Real, (Z, z) => Z * (1 - z), z => z.Re));

            Assert.That.Value(G_norm).IsEqual(0.0022682232923298762);

            var B = Enumerable.Range(0, N + 1).ToArray(i => G_norm * SpecialFunctions.BinomialCoefficient(N, i));

            Assert.That.Collection(B).ValuesAreEqual(
                0.0022682232923298762,
                0.01134111646164938,
                0.02268223292329876,
                0.02268223292329876,
                0.01134111646164938,
                0.0022682232923298762);

            // Проверяем коэффициенты передачи рассчитанного фильтра
            var H0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0).Abs;          // == (r == 1 ? 1 : Gp)
            var Hp = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fp / fd).Abs;    // == Gp
            var Hs = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fs / fd).Abs;    // <= Gs
            var Hd = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0.5).Abs;        // == 0

            Assert.That.Value(H0).IsEqual(r == 1 ? 1 : Gp, 1.59e-14);
            Assert.That.Value(Hp).IsEqual(Gp, 9.78e-15);
            Assert.That.Value(Hs).LessOrEqualsThan(Gs);
            Assert.That.Value(Hd).IsEqual(0, 5.67e-20);

            #endregion

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.I);

            Assert.That.Collection(filter.A).IsEqualTo(A);
            Assert.That.Collection(filter.B).IsEqualTo(B);
        }

        [TestMethod]
        public void TypeII_Creation()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Частота дискретизации
            const double dt = 1 / fd;           // Период дискретизации

            const double fp = 1.0 * fd / pi2;   // Граничная частота конца интервала пропускания
            const double fs = 1.6 * fd / pi2;   // Граничная частота начала интервала подавления

            const double Rp = 1.5;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 35; // Неравномерность в полосе пропускания (дБ)

            #region Аналитический расчёт

            Assert.IsTrue(fp < fs);
            Assert.IsTrue(fp < fd / 2);

            //const double wp = Consts.pi2 * fp * dt; // 0.628318530717959 рад/с
            //const double ws = Consts.pi2 * fs * dt; // 1.884955592153876 рад/с

            var eps_p = (10d.Pow(Rp / 10) - 1).Sqrt();
            var eps_s = (10d.Pow(Rs / 10) - 1).Sqrt();

            Assert.That.Value(eps_p).IsEqual(0.6422908567173865);
            Assert.That.Value(eps_s).IsEqual(56.225240418946896);

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            Assert.That.Value(Gp).IsEqual(0.84139514164519513);
            Assert.That.Value(Gs).IsEqual(0.017782794100389229);

            Assert.That.Value(Gp.In_dB()).IsEqual(-1.5, 4.67e-15);
            Assert.That.Value(Gs.In_dB()).IsEqual(-35);

            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);  // Частота подавления аналогового прототипа

            Assert.That.Value(Fp).IsEqual(869.46741682049208);
            Assert.That.Value(Fs).IsEqual(1638.7206595257194);

            var Wp = Consts.pi2 * Fp;
            var Ws = Consts.pi2 * Fs;

            var k_eps = eps_s / eps_p;
            var k_W = Ws / Wp;

            Assert.That.Value(k_eps).IsEqual(87.5385969314623);
            Assert.That.Value(k_W).IsEqual(1.8847407364824178);

            var N = (int)Math.Ceiling(arch(k_eps) / arch(k_W)); // Порядок фильтра
            Assert.That.Value(N).IsEqual(5);

            var r = N % 2;                              // Нечётность порядка фильтра
            var L = (N - r) / 2;
            var dth = Math.PI / N;                      // Угловой шаг между полюсами
            var beta = arsh(eps_s) / N;
            Assert.That.Value(beta).IsEqual(0.94451840539627485);

            var sh = Math.Sinh(beta);
            var ch = Math.Cosh(beta);
            var poles = new Complex[N];                    // Массив полюсов фильтра
            if (r != 0) poles[0] = -1 / sh;              // Если порядок фильтра нечётный, то первым добавляем центральный полюс
            for (var i = r; i < poles.Length; i += 2)   // Расчёт полюсов
            {
                var n = (i - r) / 2 + 1;
                var th = dth * (n - 0.5);

                var sin = Math.Sin(th);
                var cos = Math.Cos(th);
                var norm = (sh.Pow2() * sin.Pow2() + ch.Pow2() * cos.Pow2()).GetInverse();
                poles[i] = new Complex(-sh * sin, ch * cos) * norm;
                poles[i + 1] = poles[i].ComplexConjugate;
            }

            var zeros = new Complex[L * 2];
            for (var n = 1; n <= L; n++)
            {
                var th = dth * (n - 0.5);
                zeros[2 * n - 2] = new Complex(0, 1 / Math.Cos(th));
                zeros[2 * n - 1] = zeros[2 * n - 2].ComplexConjugate;
            }

            Assert.That.Collection(zeros).IsEqualTo(
                (0, +1.0514622242382672),
                (0, -1.0514622242382672),
                (0, +1.7013016167040798),
                (0, -1.7013016167040798)
                );

            Assert.That.Collection(poles).IsEqualTo(
                -0.916293046617597,
                (-0.16093388242115719, +0.671788117566478),
                (-0.16093388242115719, -0.671788117566478),
                (-0.5746163895448414, +0.56623918171603271),
                (-0.5746163895448414, -0.56623918171603271)
                );

            var z_zeros = r == 0
                ? DigitalFilter.ToZArray(zeros, dt, Wp)
                : DigitalFilter.ToZ(zeros, dt, Wp).Prepend(-1).ToArray();
            var z_poles = DigitalFilter.ToZArray(poles, dt, Wp);

            Assert.That.Collection(z_zeros).IsEqualTo(
                -1,
                (0.5038111428371057, +0.8638138296839021),
                (0.5038111428371057, -0.8638138296839021),
                (0.07305842913892648, +0.9973276622714083),
                (0.07305842913892648, -0.9973276622714083)
            );

            Assert.That.Collection(z_poles).IsEqualTo(
                0.33282404101389379,
                (0.65054284481357039, +0.55679574248795527),
                (0.65054284481357039, -0.55679574248795527),
                (0.44222883055670215, +0.33954724414602883),
                (0.44222883055670215, -0.33954724414602883)
            );

            var (B, b_im) = Polynom.Array.GetCoefficientsInverted(z_zeros);
            var (A, a_im) = Polynom.Array.GetCoefficientsInverted(z_poles);

            Assert.That.Collection(b_im).ElementsAreEqualTo(0, 2.23e-16);
            Assert.That.Collection(a_im).ElementsAreEqualTo(0, 1.12e-16);

            var G_norm = 1
                / (z_zeros.Aggregate(Complex.Real, (Z, z) => Z * (1 - z), z => z.Re)
                    / z_poles.Aggregate(Complex.Real, (Z, z) => Z * (1 - z), z => z.Re));

            Assert.That.Value(G_norm).IsEqual(0.033411466391200018);

            for (var i = 0; i < B.Length; i++)
                B[i] *= G_norm;

            Assert.That.Collection(B).ValuesAreEqual(
                +0.03341146639120002,
                -0.00513665024116626,
                +0.033194006484350815,
                +0.033194006484350815,
                -0.00513665024116626,
                +0.03341146639120001);

            Assert.That.Collection(A).ValuesAreEqual(
                +1,
                -2.5183673917544387,
                +2.9222427371338133,
                -1.7834584122029751,
                +0.57838133656475088,
                -0.075860624472381427);

            // Проверяем коэффициенты передачи рассчитанного фильтра
            var H0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0).Abs;          // == (r == 1 ? 1 : Gp)
            var Hp = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fp / fd).Abs;    // <= Gs  !!!
            var Hs = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fs / fd).Abs;    // <= Gs
            var Hd = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0.5).Abs;        // == 0

            Assert.That.Value(H0).IsEqual(1, 1.79e-15);
            Assert.That.Value(Hp).IsEqual(Gs, 1.29e-16);
            Assert.That.Value(Hs).LessOrEqualsThan(Gs);
            Assert.That.Value(Hd).IsEqual(0, 3.09e-18);

            #endregion

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.II);

            Assert.That.Collection(filter.A).IsEqualTo(A);
            Assert.That.Collection(filter.B).IsEqualTo(B);
        }

        [TestMethod]
        public void TypeIICorrected_Creation()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Частота дискретизации
            const double dt = 1 / fd;           // Период дискретизации

            const double fp = 1.0 * fd / pi2;   // Граничная частота конца интервала пропускания
            const double fs = 1.6 * fd / pi2;   // Граничная частота начала интервала подавления

            const double Rp = 1.5;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 35; // Неравномерность в полосе пропускания (дБ)

            #region Аналитический расчёт

            Assert.IsTrue(fp < fs);
            Assert.IsTrue(fp < fd / 2);

            //const double wp = Consts.pi2 * fp * dt; // 0.628318530717959 рад/с
            //const double ws = Consts.pi2 * fs * dt; // 1.884955592153876 рад/с

            var eps_p = (10d.Pow(Rp / 10) - 1).Sqrt();
            var eps_s = (10d.Pow(Rs / 10) - 1).Sqrt();

            Assert.That.Value(eps_p).IsEqual(0.6422908567173865);
            Assert.That.Value(eps_s).IsEqual(56.225240418946896);

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            Assert.That.Value(Gp).IsEqual(0.84139514164519513);
            Assert.That.Value(Gs).IsEqual(0.017782794100389229);

            Assert.That.Value(Gp.In_dB()).IsEqual(-1.5, 4.67e-15);
            Assert.That.Value(Gs.In_dB()).IsEqual(-35);

            var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
            var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);  // Частота подавления аналогового прототипа

            Assert.That.Value(Fp).IsEqual(869.46741682049208);
            Assert.That.Value(Fs).IsEqual(1638.7206595257194);

            var Wp = Consts.pi2 * Fp;
            var Ws = Consts.pi2 * Fs;

            var k_eps = eps_s / eps_p;
            var k_W = Ws / Wp;

            Assert.That.Value(k_eps).IsEqual(87.5385969314623);
            Assert.That.Value(k_W).IsEqual(1.8847407364824178);

            var N = (int)Math.Ceiling(arch(k_eps) / arch(k_W)); // Порядок фильтра
            Assert.That.Value(N).IsEqual(5);

            var r = N % 2;                              // Нечётность порядка фильтра
            var L = (N - r) / 2;
            var dth = Math.PI / N;                      // Угловой шаг между полюсами
            var beta = arsh(eps_s) / N;
            Assert.That.Value(beta).IsEqual(0.94451840539627485);

            var sh = Math.Sinh(beta);
            var ch = Math.Cosh(beta);
            var poles = new Complex[N];                    // Массив полюсов фильтра
            if (r != 0) poles[0] = -1 / sh;              // Если порядок фильтра нечётный, то первым добавляем центральный полюс
            for (var i = r; i < poles.Length; i += 2)   // Расчёт полюсов
            {
                var n = (i - r) / 2 + 1;
                var th = dth * (n - 0.5);

                var sin = Math.Sin(th);
                var cos = Math.Cos(th);
                var norm = (sh.Pow2() * sin.Pow2() + ch.Pow2() * cos.Pow2()).GetInverse();
                poles[i] = new Complex(-sh * sin, ch * cos) * norm;
                poles[i + 1] = poles[i].ComplexConjugate;
            }

            var zeros = new Complex[L * 2];
            for (var n = 1; n <= L; n++)
            {
                var th = dth * (n - 0.5);
                zeros[2 * n - 2] = new Complex(0, 1 / Math.Cos(th));
                zeros[2 * n - 1] = zeros[2 * n - 2].ComplexConjugate;
            }

            var kw = fp / fs;
            for (var i = 0; i < poles.Length; i++) poles[i] /= kw;
            for (var i = 0; i < zeros.Length; i++) zeros[i] /= kw;

            Assert.That.Collection(zeros).IsEqualTo(
                (0, +1.6823395587812275),
                (0, -1.6823395587812275),
                (0, +2.7220825867265277),
                (0, -2.7220825867265277)
                );

            Assert.That.Collection(poles).IsEqualTo(
                -1.4660688745881552,
                (-0.25749421187385152, +1.0748609881063649),
                (-0.25749421187385152, -1.0748609881063649),
                (-0.91938622327174624, +0.9059826907456523),
                (-0.91938622327174624, -0.9059826907456523)
                );

            var z_zeros = r == 0
                ? DigitalFilter.ToZArray(zeros, dt, Wp)
                : DigitalFilter.ToZ(zeros, dt, Wp).Prepend(-1).ToArray();
            var z_poles = DigitalFilter.ToZArray(poles, dt, Wp);

            Assert.That.Collection(z_zeros).IsEqualTo(
                -1,
                (+0.084197213369226476, +0.99644911022131843),
                (+0.084197213369226476, -0.99644911022131843),
                (-0.37722028707347804, +0.92612356358112546),
                (-0.37722028707347804, -0.92612356358112546)
            );

            Assert.That.Collection(z_poles).IsEqualTo(
                0.11054530279663551,
                (0.38604858259003627, +0.7135164941237141),
                (0.38604858259003627, -0.7135164941237141),
                (0.2009647788360637, +0.39567388309026885),
                (0.2009647788360637, -0.39567388309026885)
            );

            var (B, b_im) = Polynom.Array.GetCoefficientsInverted(z_zeros);
            var (A, a_im) = Polynom.Array.GetCoefficientsInverted(z_poles);

            Assert.That.Collection(b_im).ElementsAreEqualTo(0, 1.12e-16);
            Assert.That.Collection(a_im).ElementsAreEqualTo(0, 5.56e-17);

            var G_norm = 1
                / (z_zeros.Aggregate(Complex.Real, (Z, z) => Z * (1 - z), z => z.Re)
                    / z_poles.Aggregate(Complex.Real, (Z, z) => Z * (1 - z), z => z.Re));

            Assert.That.Value(G_norm).IsEqual(0.062095228194121721);

            for (var i = 0; i < B.Length; i++)
                B[i] *= G_norm;

            Assert.That.Collection(B).ValuesAreEqual(
                0.062095228194121721,
                0.09848589744973861,
                0.15269232505691538,
                0.15269232505691541,
                0.098485897449738652,
                0.062095228194121749);

            Assert.That.Collection(A).ValuesAreEqual(
                +1,
                -1.2845720256488355,
                +1.2951957712485003,
                -0.54541694502460425,
                +0.1756686538240842,
                -0.014328552997592952);

            // Проверяем коэффициенты передачи рассчитанного фильтра
            var H0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0).Abs;          // == (r == 1 ? 1 : Gp)
            var Hp = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fp / fd).Abs;    // >= Gp  !!!
            var Hs = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fs / fd).Abs;    // <= Gs
            var Hd = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0.5).Abs;        // == 0

            Assert.That.Value(H0).IsEqual(1, 2.23e-16);
            Assert.That.Value(Hp).GreaterOrEqualsThan(Gp);
            Assert.That.Value(Hs).LessOrEqualsThan(Gs);
            Assert.That.Value(Hd).IsEqual(0, 5.75e-18);

            #endregion

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.IICorrected);

            Assert.That.Collection(filter.A).IsEqualTo(A, 4.442e-16);
            Assert.That.Collection(filter.B).IsEqualTo(B, 1.389e-17);
        }
    }
}