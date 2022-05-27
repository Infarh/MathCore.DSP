using System;
using System.Collections;
using System.Diagnostics;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Service;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions.EllipticJacobi;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters;

[TestClassHandler("FailResultHandler")]
public class EllipticLowPass : UnitTest
{
    // ReSharper disable once UnusedMember.Local
    private static void FailResultHandler(TestResult result)
    {
        if (result.TestFailureException?.InnerException is not AssertFailedException exception) return;
        switch (exception.Data["Actual"])
        {
            case IEnumerable<Complex> actual:
                result.ToDebugEnum(actual);
                break;
            case IEnumerable actual:
                result.ToDebugEnum(actual);
                break;
            case { } actual:
                result.ToDebug(actual);
                break;
        }

        //switch (exception.Data["Expected"])
        //{
        //    case IEnumerable<Complex> expected:
        //        result.ToDebugEnum(expected);
        //        break;
        //    case IEnumerable expected:
        //        result.ToDebugEnum(expected);
        //        break;
        //    case { } expected:
        //        result.ToDebug(expected);
        //        break;
        //}
    }

    [TestMethod]
    public void Creation()
    {
        //https://ru.dsplib.org/content/filter_ellip_ap/filter_ellip_ap.html

        const double fd = 2;            // Частота дискретизации
        const double dt = 1 / fd;       // Период дискретизации

        const double pi2 = 2 * PI;

        const double fp = 2 / pi2;     // Граничная частота конца интервала пропускания
        const double fs = 4 / pi2;     // Граничная частота начала интервала подавления

        const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
        const double Rs = 40;   // Уровень подавления более 45 дБ

        var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
        var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

        #region Аналитический расчёт фильтра

        //const double wp = 2 * Math.PI * fp / fd; // 1
        //const double ws = 2 * Math.PI * fs / fd; // 1.5

        // Рассчитываем частоты цифрового фильтра
        var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);
        var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);

        (Fp, Fs).AssertEquals((0.347786966728196811, 0.9914765511533168));

        // Круговые частоты
        var Wp = Consts.pi2 * Fp;
        //var Ws = Consts.pi2 * Fs;

        // Допуск на АЧХ в интервале пропускания
        var eps_p = (Pow(10, Rp / 10) - 1).Sqrt();
        // Допуск на АЧХ в интервале подавления
        var eps_s = (Pow(10, Rs / 10) - 1).Sqrt();

        eps_p.AssertEquals(0.508847139909587520);
        eps_s.AssertEquals(99.994999874993752087);

        var kW = Fp / Fs;
        var kEps = eps_p / eps_s;

        kW.AssertEquals(0.350776794795237490);
        kEps.AssertEquals(0.005088725841749188);

        var K_w = FullEllipticIntegral(kW);
        var T_w = FullEllipticIntegralComplimentary(kW);
        var K_eps = FullEllipticIntegral(kEps);
        var T_eps = FullEllipticIntegralComplimentary(kEps);

        // Оценка снизу порядка фильтра
        var double_N = (T_eps * K_w / (K_eps * T_w));
        double_N.AssertEquals(2.776212591945528274);

        // Порядок фильтра
        int N = ((int)Ceiling(double_N)).AssertEquals(3);

        var (L, r) = N.GetDivMod(2);

        // Эллиптический модуль
        var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);

        //u.ToDebugEnum();
        u.AssertEquals(
            /*[ 0]*/ 0.3333333333333333
        );

         var m = (1 - kEps * kEps).Sqrt();
        m.AssertEquals(0.999987052350832961);

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));
        kp.AssertEquals(0.910333343015472751);

        kW = (1 - kp * kp).Sqrt();
        kW.AssertEquals(0.413875832338968130);

        var im_pz = Enumerable.Range(0, L).ToArray(i => 1 / (kW * cd_uk(u[i], kW)));

        //im_pz.ToDebugEnum();
        im_pz.AssertEquals(
            /*[ 0]*/ 2.758343343677701
        );

        var v0_complex = sn_inverse((0, 1 / eps_p), kEps) / N;
        v0_complex.AssertEquals((0, 0.30301982972239955));

        var zeros = new Complex[N - r]; // Массив нулей (на r меньше числа полюсов)
        var poles = new Complex[N];     // Массив полюсов

        // Если фильтр нечётный, то первым полюсом будет действительный полюс
        if (r != 0) poles[0] = Complex.i * sn_uk(v0_complex, kW);
        for (var i = 0; i < L; i++)
        {
            // Меняем местами действительную и мнимую часть вместо домножения на комплексную единицу
            var (p_im, p_re) = cd_uk(u[i] - v0_complex, kW);
            (poles[r + 2 * i], poles[r + 2 * i + 1]) = Complex.Conjugate(-p_re, p_im);
            (zeros[2 * i], zeros[2 * i + 1]) = Complex.Conjugate(0, 1 / (kW * cd_uk(u[i], kW)));
        }

        // Полюса
        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ -0.523721030720278202,
            /*[ 1]*/ (-0.227259770751378154, 0.976571011653282528),
            /*[ 2]*/ (-0.227259770751378154, -0.976571011653282528)
        );

        // Нули
        //zeros.ToDebugEnum();
        zeros.AssertEquals(
            /*[ 0]*/ (0, 2.758343343677700954),
            /*[ 1]*/ (0, -2.758343343677700954)
        );

        var low_pass_zeros = zeros.ToArray(p => p * Wp);
        var low_pass_poles = poles.ToArray(p => p * Wp);

        //low_pass_zeros.ToDebugEnum();
        low_pass_zeros.AssertEquals(
            /*[ 0]*/ (0, 6.027559345980697536),
            /*[ 1]*/ (0, -6.027559345980697536)
        );

        //low_pass_poles.ToDebugEnum();
        low_pass_poles.AssertEquals(
            /*[ 0]*/ -1.144440412264177143,
            /*[ 1]*/ (-0.496610314411227660, 2.134012700701830134),
            /*[ 2]*/ (-0.496610314411227660, -2.134012700701830134)
        );


        var z_zeros_enum = DigitalFilter.ToZ(low_pass_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = DigitalFilter.ToZArray(low_pass_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ -1,
            /*[ 1]*/ (-0.388513279309114723, 0.921443124560858529),
            /*[ 2]*/ (-0.388513279309114723, -0.921443124560858529)
        );

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ 0.555076812811022724,
            /*[ 1]*/ (0.452070213005796362, 0.689127155834254324),
            /*[ 2]*/ (0.452070213005796362, -0.689127155834254324)
        );

        var g_zeros = z_zeros.Multiply(z => 1 - z);
        var g_poles = z_poles.Multiply(z => 1 - z);

        var g0 = N.IsOdd() ? 1 : 1 / (1 + eps_p * eps_p).Sqrt();
        var g_norm = g0 / (g_zeros / g_poles).Abs;

        g_norm.AssertEquals(0.062093450792145871);

        var (B, im_b) = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm);
        var (A, im_a) = GetCoefficientsInverted(z_poles);

        im_b.Sum(v => v * v).AssertEquals(0, 1e-30);
        im_a.Sum(v => v * v).AssertEquals(0, 1e-30);

        // Проверяем коэффициенты передачи рассчитанного фильтра
        var H0 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0).Abs;
        var Hp = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fp / fd).Abs;
        var Hs = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fs / fd).Abs;
        var Hd = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0.5).Abs;

        const double eps = 1e-14;
        Assert.That.Value(H0).IsEqual(1, eps);                 // Коэффициент на нулевой частоте должен быть равен 1
        Assert.That.Value(Hp).GreaterOrEqualsThan(Gp, eps);    // Коэффициент передачи на граничной частоте конца интервала пропускания должен быть не меньше заданного значения Gp
        Assert.That.Value(Hs).LessOrEqualsThan(Gs);            // Коэффициент передачи на граничной частоте начала интервала подавления должен быть не больше заданного значения Gs
        Assert.That.Value(Hd).IsEqual(0, eps);                 // Коэффициент передачи на частоте, равной половине частоты дискретизации (соответствующей бесконечно большой частоте) должен быть равен нулю

        #endregion

        var filter = new DSP.Filters.EllipticLowPass(dt, fp, fs, Gp, Gs);

        Assert.That.Collection(filter.A).IsEqualTo(A, 1e-15);
        Assert.That.Collection(filter.B).IsEqualTo(B, 1e-16);
    }

    [TestMethod]
    public void Creation_fp500_fs1500_fd5000_Rp1_Rs30()
    {
        //https://ru.dsplib.org/content/filter_ellip_ap/filter_ellip_ap.html

        const double fd = 5000;         // Частота дискретизации
        const double dt = 1 / fd;       // Период дискретизации

        const double fp = 500;     // Граничная частота конца интервала пропускания
        const double fs = 1500;     // Граничная частота начала интервала подавления

        const double Rp = 1;    // Неоднородность АЧХ в интервале пропускания не более 1 дБ
        const double Rs = 30;   // Уровень подавления более 45 дБ

        var Gp = (-Rp).From_dB();   // Значения АЧХ в интервале пропускания
        var Gs = (-Rs).From_dB();   // Значения АЧХ в интервале подавления

        #region Аналитический расчёт фильтра

        // Рассчитываем частоты цифрового фильтра
        var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);
        var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);

        // Круговые частоты
        var Wp = Consts.pi2 * Fp;
        //var Ws = Consts.pi2 * Fs;

        // Допуск на АЧХ в интервале пропускания
        var eps_p = (Pow(10, Rp / 10) - 1).Sqrt();
        // Допуск на АЧХ в интервале подавления
        var eps_s = (Pow(10, Rs / 10) - 1).Sqrt();

        var k_W = Fp / Fs;
        var k_eps = eps_p / eps_s;

        var K_w = FullEllipticIntegral(k_W);
        var T_w = FullEllipticIntegralComplimentary(k_W);
        var K_eps = FullEllipticIntegral(k_eps);
        var T_eps = FullEllipticIntegralComplimentary(k_eps);

        // Оценка снизу порядка фильтра
        var double_N = T_eps * K_w / (K_eps * T_w);

        var N = (int)Ceiling(double_N); // Порядок фильтра

        var L = N / 2;  // Число комплексно сопряжённых полюсов
        var r = N % 2;  // Число (0 или 1) действительных полюсов - (чётность фильтра)

        // Эллиптический модуль
        var u = Enumerable.Range(1, L).ToArray(i => (2 * i - 1d) / N);

        var m = (1 - k_eps * k_eps).Sqrt();

        var kp = m.Power(N) * u.Aggregate(1d, (P, ui) => P * sn_uk(ui, m).Power(4));

        k_W = (1 - kp * kp).Sqrt();

        var v0_complex = sn_inverse((0, 1 / eps_p), k_eps) / N;

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

        // Рассчитываем коэффициенты полиномов числителя и знаменателя передаточной функции
        // аналогового прототипа H(p) = Q(p) / P(p)
        var analog_numerator_coefficients = GetCoefficientsInverted(zeros);     // Q(p)
        var analog_denominator_coefficients = GetCoefficientsInverted(poles);   // P(p)

        var (analog_b, _) = analog_numerator_coefficients;
        var (analog_a, _) = analog_denominator_coefficients;

        var norm_k = analog_b[^1] / analog_a[^1];

        var translated_zeros = zeros.ToArray(p => p * Wp);
        var translated_poles = poles.ToArray(p => p * Wp);

        var z_zeros = translated_zeros.ToArray(z => DigitalFilter.ToZ(z, dt));
        var z_poles = translated_poles.ToArray(z => DigitalFilter.ToZ(z, dt));

        if (r > 0)
        {
            Array.Resize(ref z_zeros, z_zeros.Length + 1);
            Array.Copy(z_zeros, 0, z_zeros, 1, z_zeros.Length - 1);
            z_zeros[0] = -1;
        }

        var k0 = (r > 0 ? 1 : 1 / (1 + eps_p * eps_p).Sqrt());

        var k_zeros = z_zeros.Multiply(z => 1 - z);
        var k_poles = z_poles.Multiply(z => 1 - z);
        //var k_zeros = z_zeros.Aggregate(Complex.Real, (Z, z) => Z * (1 - z));
        //var k_poles = z_poles.Aggregate(Complex.Real, (Z, z) => Z * (1 - z));

        var g_norm = k0 * (k_poles / k_zeros).Abs;

        var (B, Bim) = GetCoefficientsInverted(z_zeros).ToArray(b => b * g_norm);
        var (A, Aim) = GetCoefficientsInverted(z_poles);

        Bim.Sum(v => v * v).AssertEquals(0);
        Aim.Sum(v => v * v).AssertEquals(0);

        // Проверяем коэффициенты передачи рассчитанного фильтра
        var H0 = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0).Abs;
        var Hp = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fp / fd).Abs;
        var Hps = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0.5 * (fp + fs) / fd).Abs;
        var Hs = DoubleArrayDSPExtensions.FrequencyResponse(A, B, fs / fd).Abs;
        var Hd = DoubleArrayDSPExtensions.FrequencyResponse(A, B, 0.5).Abs;

        //var H = Enumerable.Range(0, 101).Select(i => fd / 2 / 100 * i)
        //   .Select(f => ($"{f,5} {DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, f / fd).Abs.In_dB().Round(2)}"))
        //   .ToArray();

        //var H0db = H0.In_dB();
        //var Hpdb = Hp.In_dB();
        //var Hpsdb = Hps.In_dB();
        //var Hsdb = Hs.In_dB();
        //var Hddb = Hd.In_dB();

        const double eps = 1e-14;
        Assert.That.Value(H0).IsEqual(Gp, eps);                // Коэффициент на нулевой частоте должен быть равен 1
        Assert.That.Value(Hp).GreaterOrEqualsThan(Gp, eps);    // Коэффициент передачи на граничной частоте конца интервала пропускания должен быть не меньше заданного значения Gp
        Assert.That.Value(Hs).LessOrEqualsThan(Gs);            // Коэффициент передачи на граничной частоте начала интервала подавления должен быть не больше заданного значения Gs
        Assert.That.Value(Hd).IsEqual(Gs, eps);                // Коэффициент передачи на частоте, равной половине частоты дискретизации (соответствующей бесконечно большой частоте) должен быть равен нулю

        #endregion

        var filter = new DSP.Filters.EllipticLowPass(dt, fp, fs, Gp, Gs);

        filter.A.AssertEquals(Accuracy.Eps(1e-15), A);
        filter.B.AssertEquals(Accuracy.Eps(1e-16), B);
    }

    [TestMethod]
    public void TransmissionCoefficient()
    {
        const double pi2 = 2 * PI;

        const double fd = 5000;      // Гц // Частота дискретизации
        const double dt = 1 / fd;    // с  // Период дискретизации
        const double fp = fd / pi2;  // Гц // Граничная частота полосы пропускания
        const double fs = 1.5 * fp;  // Гц // Граничная частота полосы запирания
        const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
        const double Rs = 45; // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.EllipticLowPass(dt, fp, fs, Gp, Gs);

        var transmission_0 = filter.GetTransmissionCoefficient(0);
        var transmission_fp = filter.GetTransmissionCoefficient(fp);
        var transmission_fs = filter.GetTransmissionCoefficient(fs);
        var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2);

        var transmission_0_abs = transmission_0.Abs;
        var transmission_fp_abs = transmission_fp.Abs;
        var transmission_fs_abs = transmission_fs.Abs;
        var transmission_fd05_abs = transmission_fd05.Abs;

        const double eps = 5.18e-14;
        Assert.That.Value(transmission_0_abs).IsEqual(Gp, eps);
        Assert.That.Value(transmission_fp_abs).IsEqual(Gp, eps);
        Assert.That.Value(transmission_fs_abs).LessOrEqualsThan(Gs);
        Assert.That.Value(transmission_fd05_abs).IsEqual(Gs, eps);
    }

    [TestMethod]
    public void ImpulseResponse()
    {
        const double pi2 = 2 * PI;

        const double fd = 5000;      // Гц // Частота дискретизации
        const double dt = 1 / fd;    // с  // Период дискретизации
        const double fp = fd / pi2;  // Гц // Граничная частота полосы пропускания
        const double fs = 1.5 * fp;  // Гц // Граничная частота полосы запирания
        const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
        const double Rs = 45; // Неравномерность в полосе пропускания (дБ)
        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.EllipticLowPass(dt, fp, fs, Gp, Gs);

        double[] expected_impulse_response =
        {
            +0.02316402765395902, +0.07890235175222579, +0.16227581269333505, +0.24853156657007530, +0.28054421394358340,
            +0.22795965026314102, +0.10322811584198478, -0.03550289833716826, -0.11558449693089311, -0.09973469008670285,
            -0.01432462412676058, +0.07033262928146783, +0.09354381942841251, +0.04804815902953240, -0.02081136468031311,
            -0.05692503754813389
        };
        var impulse_response = filter.GetImpulseResponse(expected_impulse_response.Length, 1e-10).ToArray();

        var error2 = impulse_response
           .Zip(expected_impulse_response, (a, e) => e - a)
           .Average(d => d * d);
        error2.ToDebug();

        error2.AssertEquals(0, 3.6e-3);
    }

    [TestMethod]
    public void SignalProcessing()
    {
        const double pi2 = 2 * PI;

        const double fd = 5000;      // Гц // Частота дискретизации
        const double dt = 1 / fd;    // с  // Период дискретизации
        const double fp = fd / pi2;  // Гц // Граничная частота полосы пропускания
        const double fs = 1.5 * fp;  // Гц // Граничная частота полосы запирания
        const double Rp = 1;  // Неравномерность в полосе пропускания (дБ)
        const double Rs = 45; // Неравномерность в полосе пропускания (дБ)
        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.EllipticLowPass(dt, fp, fs, Gp, Gs);

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
        Assert.That.Value(y_0.Power).IsEqual(s_0.Power, 0.21);
        // На граничной частоте сигнал должен быть ослаблен на коэффициент Gp (на Gp^2 по мощности)
        Assert.That.Value(y_fp.Power).IsEqual(s_0.Power * Gp * Gp, 2.26e-2);
        // На частоте заграждения сигнал должен быть ослаблен на коэффициент Gs (на Gs^2 по мощности)
        Assert.That.Value(y_fs.Power).LessOrEqualsThan(s_fs.Power * Gs * Gs, 2.72e-4);
        // На частоте в половину частоты дискретизации сигнал должен быть подавлен
        Assert.That.Value(y_fd05.Power).IsEqual(0, 2.26e-4);
    }
}