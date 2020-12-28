using System;
using System.Collections.Generic;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Service;

using Microsoft.VisualStudio.TestTools.UnitTesting;
// ReSharper disable InconsistentNaming
// ReSharper disable RedundantArgumentDefaultValue
// ReSharper disable ArgumentsStyleNamedExpression

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

            const double Rp = 1.5;              // Неравномерность в полосе пропускания (дБ)
            const double Rs = 35;               // Затухание в полосе подавления (дБ)

            #region Аналитический расчёт

            Assert.IsTrue(fp < fs);
            Assert.IsTrue(fp < fd / 2);

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

            var comparer = new LambdaEqualityComparer<Complex>((x1, x2) => (x1 - x2).Abs < 1e-14);

            Assert.That.Collection(poles).IsEqualTo(new Complex[]
            {
                -0.24765029577598821,
                (-0.076528150056762612, +0.979787021976867),
                (-0.076528150056762612, -0.979787021976867),
                (-0.20035329794475673, +0.605541681317744),
                (-0.20035329794475673, -0.605541681317744)
            }, comparer);

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

            const double Rp = 1.5;              // Неравномерность в полосе пропускания (дБ)
            const double Rs = 35;               // Затухание в полосе подавления (дБ)

            #region Аналитический расчёт

            Assert.IsTrue(fp < fs);
            Assert.IsTrue(fp < fd / 2);

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

            const double Rp = 1.5;              // Неравномерность в полосе пропускания (дБ)
            const double Rs = 35;               // Затухание в полосе подавления (дБ)

            #region Аналитический расчёт

            Assert.IsTrue(fp < fs);
            Assert.IsTrue(fp < fd / 2);

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

            const double kw = fp / fs;
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

            Assert.That.Value(filter.Order).IsEqual(N);
        }

        /// <summary>Тестирование коэффициента передачи фильтра Чебышева первого рода чётного порядка</summary>
        [TestMethod]
        public void TypeI_EvenOrder_TransmissionCoefficient()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.8 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 30;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.I);

            Assert.That.Value(filter.Order).IsEven();

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            var transmission_0_abs = transmission_0.Abs;
            var transmission_fp_abs = transmission_fp.Abs;
            var transmission_fs_abs = transmission_fs.Abs;
            var transmission_fd05_abs = transmission_fd05.Abs;

            const double eps = 3.42e-14;
            Assert.That.Value(transmission_0_abs).IsEqual(Gp, eps);
            Assert.That.Value(transmission_fp_abs).IsEqual(Gp, eps);
            Assert.That.Value(transmission_fs_abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05_abs).IsEqual(0, eps);
        }

        /// <summary>Тестирование коэффициента передачи фильтра Чебышева первого рода нечётного порядка</summary>
        [TestMethod]
        public void TypeI_OddOrder_TransmissionCoefficient()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.6 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 35;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.I);

            Assert.That.Value(filter.Order).IsOdd();

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            var transmission_0_abs = transmission_0.Abs;
            var transmission_fp_abs = transmission_fp.Abs;
            var transmission_fs_abs = transmission_fs.Abs;
            var transmission_fd05_abs = transmission_fd05.Abs;

            const double eps = 3.42e-14;
            Assert.That.Value(transmission_0_abs).IsEqual(1, eps);
            Assert.That.Value(transmission_fp_abs).IsEqual(Gp, eps);
            Assert.That.Value(transmission_fs_abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05_abs).IsEqual(0, eps);
        }

        [TestMethod]
        public void TypeII_EvenOrder_TransmissionCoefficient()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.8 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 30;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.II);

            Assert.That.Value(filter.Order).IsEven();

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            var transmission_0_abs = transmission_0.Abs;
            var transmission_fp_abs = transmission_fp.Abs;
            var transmission_fs_abs = transmission_fs.Abs;
            var transmission_fd05_abs = transmission_fd05.Abs;

            const double eps = 1.79e-15;
            Assert.That.Value(transmission_0_abs).IsEqual(1, eps);
            Assert.That.Value(transmission_fp_abs).LessOrEqualsThan(Gs, eps); // На граничной частоте интервала передачи фильтр имеет коэффициент передачи, равный заданному уровню затухания
            Assert.That.Value(transmission_fs_abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05_abs).LessOrEqualsThan(Gs, eps);
            // У чётного фильтра на бесконечной частоте (на fd/2)
            // коэффициент передачи равен Gs с точностью до 1e-16
        }

        [TestMethod]
        public void TypeII_OddOrder_TransmissionCoefficient()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.6 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 35;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.II);

            Assert.That.Value(filter.Order).IsOdd();

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            var transmission_0_abs = transmission_0.Abs;
            var transmission_fp_abs = transmission_fp.Abs;
            var transmission_fs_abs = transmission_fs.Abs;
            var transmission_fd05_abs = transmission_fd05.Abs;

            const double eps = 1.79e-15;
            Assert.That.Value(transmission_0_abs).IsEqual(1, eps);
            Assert.That.Value(transmission_fp_abs).LessOrEqualsThan(Gs, eps); // На граничной частоте интервала передачи фильтр имеет коэффициент передачи, равный заданному уровню затухания
            Assert.That.Value(transmission_fs_abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05_abs).IsEqual(0, eps);
            // У чётного фильтра на бесконечной частоте (на fd/2)
            // коэффициент передачи равен 0 с точностью до 1e-16
        }

        [TestMethod]
        public void TypeIICorrected_EvenOrder_TransmissionCoefficient()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.8 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 30;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.IICorrected);

            Assert.That.Value(filter.Order).IsEven();

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            var transmission_0_abs = transmission_0.Abs;
            var transmission_fp_abs = transmission_fp.Abs;
            var transmission_fs_abs = transmission_fs.Abs;
            var transmission_fd05_abs = transmission_fd05.Abs;

            const double eps = 3.42e-14;
            Assert.That.Value(transmission_0_abs).IsEqual(1, eps);
            Assert.That.Value(transmission_fp_abs).GreaterOrEqualsThan(Gp);
            Assert.That.Value(transmission_fs_abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05_abs).LessOrEqualsThan(Gs, eps);
        }

        [TestMethod]
        public void TypeIICorrected_OddOrder_TransmissionCoefficient()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.6 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 35;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.IICorrected);

            Assert.That.Value(filter.Order).IsOdd();

            var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
            var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
            var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
            var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

            var transmission_0_abs = transmission_0.Abs;
            var transmission_fp_abs = transmission_fp.Abs;
            var transmission_fs_abs = transmission_fs.Abs;
            var transmission_fd05_abs = transmission_fd05.Abs;

            const double eps = 3.42e-14;
            Assert.That.Value(transmission_0_abs).IsEqual(1, eps);
            Assert.That.Value(transmission_fp_abs).GreaterOrEqualsThan(Gp);
            Assert.That.Value(transmission_fs_abs).LessOrEqualsThan(Gs);
            Assert.That.Value(transmission_fd05_abs).IsEqual(0, eps);
        }

        /// <summary>Тестирование импульсной характеристики фильтра Чебышева первого рода чётного порядка</summary>
        [TestMethod]
        public void TypeI_EvenOrder_ImpulseResponse()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.8 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 30;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.I);

            Assert.That.Value(filter.Order).IsEven();

            double[] expected_impulse_response =
            {
                +0.00884055732125126,   //  0
                +0.05614365197974830,   //  1
                +0.16099203897373400,   //  2
                +0.27575709786712000,   //  3
                +0.30792640296908800,   //  4
                +0.21418298648174700,   //  5
                +0.04842642349831240,   //  6
                -0.08470720459467050,   //  7
                -0.11573677212917000,   //  8
                -0.05772015042422970,   //  9
                +0.01790936028253510,   // 10
                +0.04635399561641480,   // 11
                +0.01696035720144810,   // 12
                -0.03089992850593000,   // 13
                -0.05045951649249290,   // 14
                -0.02740883863154980    // 15
            };

            var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();

            const double eps = 9.4e-7;
            Assert.That.Collection(impulse_response)
               .IsEqualTo(expected_impulse_response, eps);
        }

        /// <summary>Тестирование импульсной характеристики фильтра Чебышева первого рода нечётного порядка</summary>
        [TestMethod]
        public void TypeI_OddOrder_ImpulseResponse()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.6 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 35;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.I);

            Assert.That.Value(filter.Order).IsOdd();

            double[] expected_impulse_response =
            {
                +0.00226822329232988,   //  0
                +0.01839598031217860,   //  1
                +0.06923715452579570,   //  2
                +0.16078383912296200,   //  3
                +0.25625891861547400,   //  4
                +0.28937995857462800,   //  5
                +0.21928291005052700,   //  6
                +0.07371124356918550,   //  7
                -0.06457862559364690,   //  8
                -0.11880487485266500,   //  9
                -0.07442276580731260,   // 10
                +0.01561090045052300,   // 11
                +0.07749926275065900,   // 12
                +0.07349323870675280,   // 13
                +0.02314957138213890,   // 14
                -0.02296390124555760,   // 15
            };
            // impulse_response.Select(s => s.ToString("+0.00000000000000000',';-0.00000000000000000','", System.Globalization.CultureInfo.InvariantCulture)).Select((s,i) => $"{s}   // {i,2}").ToSeparatedStr("\r\n")
            var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();

            const double eps = 9.4e-7;
            Assert.That.Collection(impulse_response)
               .IsEqualTo(expected_impulse_response, eps);
        }

        /// <summary>Тестирование импульсной характеристики фильтра Чебышева второго рода чётного порядка</summary>
        [TestMethod]
        public void TypeII_EvenOrder_ImpulseResponse()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.8 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 30;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.II);

            Assert.That.Value(filter.Order).IsEven();

            double[] expected_impulse_response =
            {
                +0.05012504332253210,   //  0
                +0.09392699112998970,   //  1
                +0.16880264674281600,   //  2
                +0.20428669020597100,   //  3
                +0.22480768572184000,   //  4
                +0.20594716920198400,   //  5
                +0.14351733262540600,   //  6
                +0.05938964294902360,   //  7
                -0.01582220884364540,   //  8
                -0.06007473056575070,   //  9
                -0.06779542575678730,   // 10
                -0.04802429119514430,   // 11
                -0.01693649613045550,   // 12
                +0.01029315475076850,   // 13
                +0.02479638209856960,   // 14
                +0.02529890931233970,   // 15
            };

            var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();

            const double eps = 9.4e-7;
            Assert.That.Collection(impulse_response)
               .IsEqualTo(expected_impulse_response, eps);
        }

        /// <summary>Тестирование импульсной характеристики фильтра Чебышева второго рода нечётного порядка</summary>
        [TestMethod]
        public void TypeII_OddOrder_ImpulseResponse()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.6 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 35;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.II);

            Assert.That.Value(filter.Order).IsOdd();

            double[] expected_impulse_response =
            {
                +0.03341146639120000,   //  0
                +0.07900569722913120,   //  1
                +0.13452296315034400,   //  2
                +0.20068638620376200,   //  3
                +0.22873545544693900,   //  4
                +0.21975234841160600,   //  5
                +0.17110029908928400,   //  6
                +0.09079561772546840,   //  7
                +0.00350716717435198,   //  8
                -0.06109287376964880,   //  9
                -0.08446355604367730,   // 10
                -0.06746188706588570,   // 11
                -0.02716807079979880,   // 12
                +0.01368462023930950,   // 13
                +0.03775672984570410,   // 14
                +0.03925364904843500,   // 15
            };
            // impulse_response.Select(s => s.ToString("+0.00000000000000000',';-0.00000000000000000','", System.Globalization.CultureInfo.InvariantCulture)).Select((s,i) => $"{s}   // {i,2}").ToSeparatedStr("\r\n")
            var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();

            const double eps = 9.4e-7;
            Assert.That.Collection(impulse_response)
               .IsEqualTo(expected_impulse_response, eps);
        }

        /// <summary>Тестирование импульсной характеристики фильтра Чебышева второго рода чётного порядка</summary>
        [TestMethod]
        public void TypeIICorrected_EvenOrder_ImpulseResponse()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.8 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 30;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.IICorrected);

            Assert.That.Value(filter.Order).IsEven();

            double[] expected_impulse_response =
            {
                +0.09748158833651720,   //  0
                +0.24801587857989800,   //  1
                +0.36603707048083000,   //  2
                +0.31591693100205300,   //  3
                +0.13960087531424400,   //  4
                -0.06424082484097170,   //  5
                -0.12598772097936600,   //  6
                -0.04986798951161600,   //  7
                +0.03782993452902130,   //  8
                +0.05321816100665740,   //  9
                +0.01375504036863090,   // 10
                -0.02100526869452690,   // 11
                -0.02161750176356920,   // 12
                -0.00237250156111623,   // 13
                +0.01065544118382300,   // 14
                +0.00836193578789439,   // 15
            };

            var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();

            const double eps = 9.4e-7;
            Assert.That.Collection(impulse_response)
               .IsEqualTo(expected_impulse_response, eps);
        }

        /// <summary>Тестирование импульсной характеристики фильтра Чебышева второго рода нечётного порядка</summary>
        [TestMethod]
        public void TypeIICorrected_OddOrder_ImpulseResponse()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.6 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 35;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.IICorrected);

            Assert.That.Value(filter.Order).IsOdd();

            double[] expected_impulse_response =
            {
                +0.06209522819412170,   //  0
                +0.17825169051418800,   //  1
                +0.30124398324431900,   //  2
                +0.34265887271796400,   //  3
                +0.23479927381146400,   //  4
                +0.05378155735164640,   //  5
                -0.09849782313435830,   //  6
                -0.12399972373459100,   //  7
                -0.03871630611893290,   //  8
                +0.05106425231274820,   //  9
                +0.06618294738144070,   // 10
                +0.01813356298180970,   // 11
                -0.02955019258530920,   // 12
                -0.03487370097223260,   // 13
                -0.00752873612681361,   // 14
                +0.01714269771629890,   // 15
            };
            // impulse_response.Select(s => s.ToString("+0.00000000000000000',';-0.00000000000000000','", System.Globalization.CultureInfo.InvariantCulture)).Select((s,i) => $"{s}   // {i,2}").ToSeparatedStr("\r\n")
            var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();

            const double eps = 9.4e-7;
            Assert.That.Collection(impulse_response)
               .IsEqualTo(expected_impulse_response, eps);
        }

        [TestMethod]
        public void TypeI_EvenOrder_SignalProcessing()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.8 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 30;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.I);

            Assert.That.Value(filter.Order).IsEven();

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
            Assert.That.Value(y_0.Power).IsEqual(s_0.Power * Gp * Gp, 1.45e-3);
            // На граничной частоте сигнал должен быть ослаблен на коэффициент Gp (на Gp^2 по мощности)
            Assert.That.Value(y_fp.Power).IsEqual(s_fp.Power * Gp * Gp, 8.24e-3);
            // На частоте заграждения сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
            Assert.That.Value(y_fs.Power).LessOrEqualsThan(s_fs.Power * Gs * Gs);
            // На частоте в половину частоты дискретизации сигнал должен быть подавлен
            Assert.That.Value(y_fd05.Power).IsEqual(0, 1.55e-4);
        }

        [TestMethod]
        public void TypeI_OddOrder_SignalProcessing()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.5 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 35;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.I);

            Assert.That.Value(filter.Order).IsOdd();

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
            Assert.That.Value(y_0.Power).IsEqual(1, 5.65e-3);
            // На граничной частоте сигнал должен быть ослаблен на коэффициент Gp (на Gp^2 по мощности)
            Assert.That.Value(y_fp.Power).IsEqual(s_fp.Power * Gp * Gp, 1.34e-2);
            // На частоте заграждения сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
            Assert.That.Value(y_fs.Power).LessOrEqualsThan(s_fs.Power * Gs * Gs, 5.61e-5);
            // На частоте в половину частоты дискретизации сигнал должен быть подавлен
            Assert.That.Value(y_fd05.Power).IsEqual(0, 1.55e-4);
        }

        [TestMethod]
        public void TypeII_EvenOrder_SignalProcessing()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.8 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 30;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.II);

            Assert.That.Value(filter.Order).IsEven();

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
            Assert.That.Value(y_0.Power).IsEqual(1, 3.35e-3);
            // На граничной частоте сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
            Assert.That.Value(y_fp.Power).IsEqual(s_fp.Power * Gs * Gs, 8.24e-3);
            // На частоте заграждения сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
            Assert.That.Value(y_fs.Power).LessOrEqualsThan(s_fs.Power * Gs * Gs);
            // На частоте в половину частоты дискретизации сигнал должен быть подавлен
            Assert.That.Value(y_fd05.Power).IsEqual(y_fd05.Power * Gs * Gs, 2.11e-3);
        }

        [TestMethod]
        public void TypeII_OddOrder_SignalProcessing()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.5 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 35;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.II);

            Assert.That.Value(filter.Order).IsOdd();

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
            Assert.That.Value(y_0.Power).IsEqual(1, 3.59e-3);
            // На граничной частоте сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
            Assert.That.Value(y_fp.Power).IsEqual(s_fp.Power * Gs * Gs, 3.56e-4);
            // На частоте заграждения сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
            Assert.That.Value(y_fs.Power).LessOrEqualsThan(s_fs.Power * Gs * Gs, 1.77e-4);
            // На частоте в половину частоты дискретизации сигнал должен быть подавлен
            Assert.That.Value(y_fd05.Power).IsEqual(0, 1.16e-4);
        }

        [TestMethod]
        public void TypeIICorrected_EvenOrder_SignalProcessing()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.8 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 30;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.IICorrected);

            Assert.That.Value(filter.Order).IsEven();

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
            Assert.That.Value(y_0.Power).IsEqual(1, 1.87e-3);
            // На граничной частоте сигнал должен быть ослаблен на коэффициент Gp (на Gp^2 по мощности)
            Assert.That.Value(y_fp.Power).GreaterOrEqualsThan(s_fp.Power * Gp * Gp, 6.57e-2);
            // На частоте заграждения сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
            Assert.That.Value(y_fs.Power).LessOrEqualsThan(s_fs.Power * Gs * Gs, 7.28e-5);
            // На частоте в половину частоты дискретизации сигнал должен быть подавлен
            Assert.That.Value(y_fd05.Power).IsEqual(y_fd05.Power * Gs * Gs, 2.20e-3);
        }

        [TestMethod]
        public void TypeIICorrected_OddOrder_SignalProcessing()
        {
            const double pi2 = 2 * Math.PI;

            const double fd = 5000;             // Гц // Частота дискретизации
            const double dt = 1 / fd;           // с  // Период дискретизации
            const double fp = 1.0 * fd / pi2;   // Гц // Граничная частота полосы пропускания
            const double fs = 1.5 * fd / pi2;   // Гц // Граничная частота полосы запирания
            const double Rp = 1.5;              // дБ // Неравномерность в полосе пропускания
            const double Rs = 35;               // дБ // Затухание в полосе подавления

            var Gp = (-Rp).From_dB();
            var Gs = (-Rs).From_dB();

            var filter = new ChebyshevLowPass(fp, fs, dt, Gp, Gs, ChebyshevLowPass.ChebyshevType.IICorrected);

            Assert.That.Value(filter.Order).IsOdd();

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
            Assert.That.Value(y_0.Power).IsEqual(1, 2.40e-3);
            // На граничной частоте сигнал должен быть ослаблен на коэффициент Gp (на Gp^2 по мощности)
            Assert.That.Value(y_fp.Power).GreaterOrEqualsThan(s_fp.Power * Gp * Gp, 1.68e-1);
            // На частоте заграждения сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
            Assert.That.Value(y_fs.Power).LessOrEqualsThan(s_fs.Power * Gs * Gs, 1.67e-4);
            // На частоте в половину частоты дискретизации сигнал должен быть подавлен
            Assert.That.Value(y_fd05.Power).IsEqual(0, 1.74e-4);
        }
    }
}