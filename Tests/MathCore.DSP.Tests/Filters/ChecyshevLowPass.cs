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

            Assert.IsTrue(fp < fs);
            Assert.IsTrue(fp < fd / 2);

            //const double wp = Consts.pi2 * fp * dt; // 0.628318530717959 рад/с
            //const double ws = Consts.pi2 * fs * dt; // 1.884955592153876 рад/с

            const double Rp = 1.5;  // Неравномерность в полосе пропускания (дБ)
            const double Rs = 35; // Неравномерность в полосе пропускания (дБ)

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
        }
    }
}