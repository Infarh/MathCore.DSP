using System;
using System.Linq;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Linq.Enumerable;
using static System.Math;
using static MathCore.SpecialFunctions.EllipticJacobi;
// ReSharper disable ArgumentsStyleLiteral

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class EllipticHighPass : UnitTest
{
    [TestMethod]
    public void Creation()
    {
        const double fd = 5000;                 // Гц // Частота дискретизации
        const double fs = 0.5 * fd / PI;   // Гц // Граничная частота полосы пропускания
        const double fp = 1.0 * fd / PI;   // Гц // Граничная частота полосы запирания
        const double dt = 1 / fd;               // 2с // Период дискретизации
        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 30;                   // Неравномерность в полосе пропускания (дБ)
        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        const double ws = 2 * PI * fs / fd;
        const double wp = 2 * PI * fp / fd;

        var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);
        var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);

        var Ws = 2 * PI * Fs;
        var Wp = 2 * PI * Fp;

        var eps_p = Sqrt(Pow(10, Rp / 10) - 1);
        var eps_s = Sqrt(Pow(10, Rs / 10) - 1);

        var k_eps = eps_p / eps_s;
        k_eps.AssertEquals(0.016099211048699186);
        //var k_W = Ws / Wp;

        var k_w = Min(wp, ws) / Max(wp, ws);
        k_w.AssertEquals(0.5);

        var Kw = FullEllipticIntegral(k_w);
        var Tw = FullEllipticIntegralComplimentary(k_w);
        var K_eps = FullEllipticIntegral(k_eps);
        var T_eps = FullEllipticIntegralComplimentary(k_eps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N);
        N.AssertEquals(3);

        var L = N / 2; // Число парных (комплексно-сопряжённых) полюсов
        var r = N % 2; // Число нечётных полюсов (1, или 0)

        L.AssertEquals(1);
        r.AssertEquals(1);

        var u = Range(1, L).ToArray(i => (2 * i - 1d) / N);

        var m = (1 - k_eps * k_eps).Sqrt();
        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));

        k_w = (1 - kp * kp).Sqrt();

        //var im_pz = Range(0, L).ToArray(i => 1 / (k_w * cd_uk(u[i], k_w)));

        var v0_complex = sn_inverse((0, 1 / eps_p), k_eps) / N;

        var zeros = new Complex[N - r]; // Массив нулей (на r меньше числа полюсов - либо на 1, либо на 0)
        var poles = new Complex[N];     // Массив полюсов

        // Если фильтр нечётный, то первым полюсом будет действительный полюс
        if (r != 0) poles[0] = Complex.i * sn_uk(v0_complex, k_w);
        for (var i = 0; i < L; i++)
        {
            // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, k_w);

            poles[r + 2 * i] = new(-p_re, p_im);
            poles[r + 2 * i + 1] = new(-p_re, -p_im);

            var p0_im = 1 / (k_w * cd_uk(u[i], k_w));
            zeros[2 * i] = Complex.ImValue(p0_im);
            zeros[2 * i + 1] = Complex.ImValue(-p0_im);
        }

        zeros.AssertEquals(
            (0, 1.9535902076550578),
            (0, -1.9535902076550578));

        poles.AssertEquals(
            -0.5595579001364851,
            (-0.20528349152310388, +0.98694709202692321),
            (-0.20528349152310388, -0.98694709202692321));

        var translated_zeros = AnalogBasedFilter.TransformToHighPass(zeros, Fp);
        var translated_poles = AnalogBasedFilter.TransformToHighPass(poles, Fp);

        if (r > 0)
            translated_zeros = translated_zeros.Prepend(0);

        translated_zeros.AssertEquals(
            0,
            (0, -7972.0287220537257),
            (0, +7972.0287220537257));

        translated_poles.AssertEquals(
            -27832.825240694943,
            (-3146.1154034951815, -15125.665613063766),
            (-3146.1154034951815, +15125.665613063766));

        var z_zeros = DigitalFilter.ToZArray(translated_zeros, dt);
        var z_poles = DigitalFilter.ToZArray(translated_poles, dt);

        z_zeros.AssertEquals(
            1,
            (0.2228433849507625, -0.97485425874008791),
            (0.2228433849507625, +0.97485425874008791));

        z_poles.AssertEquals(
            -0.471358539237853,
            (-0.34532293898194616, -0.753258738841563),
            (-0.34532293898194616, +0.753258738841563));

        var Gz0 = (z_zeros.Multiply(z => 1 + z) / z_poles.Multiply(z => 1 + z)).Abs;
        Gz0.AssertEquals(9.2898762575228027);

        var G_norm = (r > 0 ? 1 : Gp) / Gz0;
        G_norm.AssertEquals(G_norm);

        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToArray(v => (v * G_norm).Re);
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        B.AssertEquals(
            +0.1076440602952289,
            -0.15561959382729446,
            +0.15561959382729446,
            -0.1076440602952289);

        A.AssertEquals(
            1,
            1.1620044172017452,
            1.0121884919960158,
            0.3236567665492236);

        var H0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0, dt).Abs;
        var Hs = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fs, dt).Abs;
        var Hp = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fp, dt).Abs;
        var Hd = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fd / 2, dt).Abs;

        H0.AssertEquals(0);
        Hs.AssertEquals(Gs, Eps: 9.48e-4);
        Hp.AssertEquals(Gp, Eps: 1e-15);
        Hd.AssertEquals(1, Eps: 1e-15);

        /* ------------------------------------------------------------------------- */

        var filter = new DSP.Filters.EllipticHighPass(dt, fs, fp, Gp, Gs);

        var actual_B = filter.B;
        var actual_A = filter.A;

        var comparer = GetComparer();
        actual_B.AssertThatCollection().IsEqualTo(B, comparer);
        actual_A.AssertThatCollection().IsEqualTo(A, comparer);

        //actual_B.AssertEquals(B);
        //actual_A.AssertEquals(A);

        //var transmission_0 = filter.GetTransmissionCoefficient(0, dt);
        //var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
        //var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
        //var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

        //var transmission_0_abs = transmission_0.Abs;
        //var transmission_fp_abs = transmission_fp.Abs;
        //var transmission_fs_abs = transmission_fs.Abs;
        //var transmission_fd05_abs = transmission_fd05.Abs;

        //Assert.That.Value(transmission_0_abs).IsEqual(0, 1e-15);
        //Assert.That.Value(transmission_fp_abs).LessOrEqualsThan(Gs, 1e-15);
        //Assert.That.Value(transmission_fs_abs).IsEqual(Gp);
        //Assert.That.Value(transmission_fd05_abs).IsEqual(1, 1e-15);
    }
}