using System;
using System.Diagnostics;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Linq.Enumerable;
using static System.Math;

using static MathCore.DSP.DoubleArrayDSPExtensions;
using static MathCore.DSP.Filters.DigitalFilter;
using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions;

// ReSharper disable InconsistentNaming
// ReSharper disable RedundantArgumentDefaultValue
// ReSharper disable ArgumentsStyleNamedExpression

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ChebyshevHighPass : ChebyshevFiltersTests
{
    [TestMethod]
    public void TypeI_Even_Creation()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        const double ws = Consts.pi2 * fs;
        const double wp = Consts.pi2 * fp;

        var Fp = ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
        //var Fs = ToAnalogFrequency(fs, dt);  // Частота подавления аналогового прототипа

        //(Fp, Fs).ToDebug();

        var Wp = Consts.pi2 * Fp;
        //var Ws = Consts.pi2 * Fs;

        var eps_p = Sqrt(Pow(10, Rp / 10) - 1);
        var eps_s = Sqrt(Pow(10, Rs / 10) - 1);

        var kEps = eps_s / eps_p;

        const double kW = wp / ws;
        var k_eps = eps_s / eps_p;

        k_eps.AssertEquals(196.512846456719728394);
        kW.AssertEquals(2);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(5);

        var beta = arcsh(1 / eps_p) / N;
        beta.AssertEquals(0.285595071772705045);

        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -Sinh(beta);
        var r = N % 2;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (im, re) = Complex.Trigonometry.Cos(new(dth * (i - r + 1), -beta));
            (poles[i], poles[i + 1]) = Complex.Conjugate(-re, im);
        }

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/  -0.289493341235612878,
            /*[ 1]*/ (-0.089458362200190114, +0.990107112003389300),
            /*[ 2]*/ (-0.089458362200190114, -0.990107112003389300),
            /*[ 3]*/ (-0.234205032817996567, +0.611919847721093646),
            /*[ 4]*/ (-0.234205032817996567, -0.611919847721093646)
        );

        //Wp.ToDebug();
        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/  -14.004469646415172335,
            /*[ 1]*/ (-00.366970242550695547, -4.061552638645113511),
            /*[ 2]*/ (-00.366970242550695547, +4.061552638645113511),
            /*[ 3]*/ (-02.211783975780333922, -5.778844704431379320),
            /*[ 4]*/ (-02.211783975780333922, +5.778844704431379320)
        );

        //var z_zeros = ToZArray(high_pass_zeros, dt);
        var z_poles = ToZArray(high_pass_poles, dt);

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/  0.176315949518621284,
            /*[ 1]*/ (0.888848859250257473, -0.376671590174077864),
            /*[ 2]*/ (0.888848859250257473, +0.376671590174077864),
            /*[ 3]*/ (0.686677174937777801, -0.438823170219132586),
            /*[ 4]*/ (0.686677174937777801, +0.438823170219132586)
        );

        //var k_zeros = z_zeros.Multiply(z => 1 + z).Re;
        var k_poles = z_poles.Multiply(z => (1 + z) / 2).Re;

        var g_norm = N.IsEven()
            ? Gp * k_poles
            : 1 * k_poles;

        //g_norm.ToDebug();

        //var zz0 = Complex.Exp(-Consts.pi2 * 0 / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fs / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fp / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fd / 2 / fd);
        //var P0 = (1 - zz0).Power(N);
        //var Pp = z_poles.Multiply(z => 1 - z * zz0);
        //var k0 = g_norm * P0 / Pp;
        //k0.Abs.ToDebug();

        //var B = GetCoefficientsInverted(z_zeros).ToArray(z => z.Re * g_norm);

        var B = new double[N + 1];
        for (var i = 0; i < B.Length; i++)
            B[i] = BinomialCoefficient(N, i) * (i % 2 == 0 ? g_norm : -g_norm);

        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        B.AssertEquals(
            /*[ 0]*/ +0.4142030534319225,
            /*[ 1]*/ -2.0710152671596123,
            /*[ 2]*/ +4.1420305343192245,
            /*[ 3]*/ -4.1420305343192245,
            /*[ 4]*/ +2.0710152671596123,
            /*[ 5]*/ -0.4142030534319225
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ +1,
            /*[ 1]*/ -3.32736801789469200,
            /*[ 2]*/ +4.59301473062132500,
            /*[ 3]*/ -3.17229294030296400,
            /*[ 4]*/ +1.05270199595006000,
            /*[ 5]*/ -0.10912002505247781
        );

        //Debug.WriteLine("---------------------------------------------");

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

        filter.B.AssertEquals(Accuracy.Eps(1e-14), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-14), A);

        //var h_0 = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, 0, dt);
        //var h_p = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fs, dt);
        //var h_s = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fp, dt);
        //var h_d = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fd / 2, dt);

        //filter.GetTransmissionCoefficient(0).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fs).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fp).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fd / 2).Abs.ToDebug();

        //var h_0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0, dt);
        //var h_p = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fs, dt);
        //var h_s = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fp, dt);
        //var h_d = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fd / 2, dt);

        //h_0.Abs.ToDebug();
        //h_p.Abs.ToDebug();
        //h_s.Abs.ToDebug();
        //h_d.Abs.ToDebug();
    }

    [TestMethod]
    public void TypeI_Odd_Creation()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        const double ws = Consts.pi2 * fs;
        const double wp = Consts.pi2 * fp;

        var Fp = ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
        //var Fs = ToAnalogFrequency(fs, dt);  // Частота подавления аналогового прототипа

        //(Fp, Fs).ToDebug();

        var Wp = Consts.pi2 * Fp;
        //var Ws = Consts.pi2 * Fs;

        var eps_p = Sqrt(Pow(10, Rp / 10) - 1);
        var eps_s = Sqrt(Pow(10, Rs / 10) - 1);

        var kEps = eps_s / eps_p;

        const double kW = wp / ws;
        var k_eps = eps_s / eps_p;

        k_eps.AssertEquals(196.512846456719728394);
        kW.AssertEquals(1.333333333333333259);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(8);

        var beta = arcsh(1 / eps_p) / N;
        beta.AssertEquals(0.178496919857940667);

        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -Sinh(beta);
        var r = N % 2;
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (im, re) = Complex.Trigonometry.Cos(new(dth * (i - r + 1), -beta));
            (poles[i], poles[i + 1]) = Complex.Conjugate(-re, im);
        }

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.035008233302118830, +0.996451282625766654),
            /*[ 1]*/ (-0.035008233302118830, -0.996451282625766654),
            /*[ 2]*/ (-0.099695013736534027, +0.844750607699364164),
            /*[ 3]*/ (-0.099695013736534027, -0.844750607699364164),
            /*[ 4]*/ (-0.149204132067111944, +0.564444310434061181),
            /*[ 5]*/ (-0.149204132067111944, -0.564444310434061181),
            /*[ 6]*/ (-0.175998273829297280, +0.198206483605588013),
            /*[ 7]*/ (-0.175998273829297280, -0.198206483605588013)
        );

        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/ (-00.142766913488235864, -04.063623343519719988),
            /*[ 1]*/ (-00.142766913488235864, +04.063623343519719988),
            /*[ 2]*/ (-00.558617485721441875, -04.733360705298244930),
            /*[ 3]*/ (-00.558617485721441875, +04.733360705298244930),
            /*[ 4]*/ (-01.774643033921778246, -06.713535005170165881),
            /*[ 5]*/ (-01.774643033921778246, +06.713535005170165881),
            /*[ 6]*/ (-10.155424345333114999, -11.436878926229493203),
            /*[ 7]*/ (-10.155424345333114999, +11.436878926229493203)
        );

        var z_poles = ToZArray(high_pass_poles, dt);

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.908163284541766314, -0.384954902154918066),
            /*[ 1]*/ (0.908163284541766314, +0.384954902154918066),
            /*[ 2]*/ (0.847710439616761957, -0.425411873912513561),
            /*[ 3]*/ (0.847710439616761957, +0.425411873912513561),
            /*[ 4]*/ (0.677531966042908551, -0.517214888840021314),
            /*[ 5]*/ (0.677531966042908551, +0.517214888840021314),
            /*[ 6]*/ (0.159654729315091504, -0.439815755318950008),
            /*[ 7]*/ (0.159654729315091504, +0.439815755318950008)
        );

        var k_poles = z_poles.Multiply(z => (1 + z) / 2).Re;

        var g_norm = N.IsEven()
            ? Gp * k_poles
            : 1 * k_poles;

        //g_norm.ToDebug();

        //var zz0 = Complex.Exp(-Consts.pi2 * 0 / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fs / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fp / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fd / 2 / fd);
        //var P0 = (1 - zz0).Power(N);
        //var Pp = z_poles.Multiply(z => 1 - z * zz0);
        //var k0 = g_norm * P0 / Pp;
        //k0.Abs.ToDebug();

        var B = new double[N + 1];
        for (var i = 0; i < B.Length; i++)
            B[i] = BinomialCoefficient(N, i) * (i % 2 == 0 ? g_norm : -g_norm);

        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        B.AssertEquals(
            /*[ 0]*/ +00.22481197152185653,
            /*[ 1]*/ -01.79849577217485220,
            /*[ 2]*/ +06.29473520261198250,
            /*[ 3]*/ -12.58947040522396500,
            /*[ 4]*/ +15.73683800652995800,
            /*[ 5]*/ -12.58947040522396500,
            /*[ 6]*/ +06.29473520261198250,
            /*[ 7]*/ -01.79849577217485220,
            /*[ 8]*/ +00.22481197152185653
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ +01,
            /*[ 1]*/ -05.18612083903305700,
            /*[ 2]*/ +12.21012661912232400,
            /*[ 3]*/ -16.94342012085136400,
            /*[ 4]*/ +15.21333307296616700,
            /*[ 5]*/ -09.16724322052108300,
            /*[ 6]*/ +03.72978925944532720,
            /*[ 7]*/ -00.98499942579790050,
            /*[ 8]*/ +00.13922172591859294
        );

        //Debug.WriteLine("---------------------------------------------");

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

        filter.B.AssertEquals(Accuracy.Eps(1e-14), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-14), A);

        //var h_0 = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, 0, dt);
        //var h_p = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fs, dt);
        //var h_s = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fp, dt);
        //var h_d = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fd / 2, dt);

        //filter.GetTransmissionCoefficient(0).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fs).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fp).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fd / 2).Abs.ToDebug();

        //var h_0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0, dt);
        //var h_p = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fs, dt);
        //var h_s = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fp, dt);
        //var h_d = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fd / 2, dt);

        //h_0.Abs.ToDebug();
        //h_p.Abs.ToDebug();
        //h_s.Abs.ToDebug();
        //h_d.Abs.ToDebug();
    }

    [TestMethod]
    public void TypeII_Even_Creation()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        const double ws = Consts.pi2 * fs;
        const double wp = Consts.pi2 * fp;

        var Fp = ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
        //var Fs = ToAnalogFrequency(fs, dt);  // Частота подавления аналогового прототипа

        //(Fp, Fs).ToDebug();

        var Wp = Consts.pi2 * Fp;
        //var Ws = Consts.pi2 * Fs;

        var eps_p = Sqrt(Pow(10, Rp / 10) - 1);
        var eps_s = Sqrt(Pow(10, Rs / 10) - 1);

        var kEps = eps_s / eps_p;

        const double kW = wp / ws;
        var k_eps = eps_s / eps_p;

        k_eps.AssertEquals(196.512846456719728394);
        kW.AssertEquals(2);

        var N = (int)Ceiling(arcch(kEps) / arcch(kW)); // Порядок фильтра
        N.AssertEquals(5);

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(1.059658473122097044);

        var r = N % 2;
        var poles = new Complex[N];
        var zeros = new Complex[N - r];
        if (N.IsOdd())
            poles[0] = -1 / Sinh(beta);
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var th = dth * (i - r + 1);
            var (re, im) = 1 / Complex.Trigonometry.Cos(new(th, beta));
            (poles[i], poles[i + 1]) = Complex.Conjugate(-im, re);

            (zeros[i - r], zeros[i]) = Complex.Conjugate(0, 1 / Cos(th));
        }

        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ +1.0514622242382672,
            /*[ 1]*/ -1.0514622242382672,
            /*[ 2]*/ +1.7013016167040798,
            /*[ 3]*/ -1.7013016167040798
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/  -0.787770266856923418,
            /*[ 1]*/ (-0.155915595278454944, +0.610870317639875871),
            /*[ 2]*/ (-0.155915595278454944, -0.610870317639875871),
            /*[ 3]*/ (-0.524799478613158787, +0.485389011298822448),
            /*[ 4]*/ (-0.524799478613158787, -0.485389011298822448)
        );

        var high_pass_zeros = AnalogBasedFilter.TransformToHighPassW(zeros, Wp);
        if (N.IsOdd())
            high_pass_zeros = high_pass_zeros.Prepend(0);
        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        //high_pass_zeros.ToIm().ToDebugEnum();
        high_pass_zeros.ToIm().AssertEquals(
            /*[ 0]*/ +0,
            /*[ 1]*/ -3.8557740037789000,
            /*[ 2]*/ +3.8557740037789000,
            /*[ 3]*/ -2.3829993872736255,
            /*[ 4]*/ +2.3829993872736255
        );

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/  -5.146425145428575654,
            /*[ 1]*/ (-1.590332830215989812, -6.230852785522531079),
            /*[ 2]*/ (-1.590332830215989812, +6.230852785522531079),
            /*[ 3]*/ (-4.163545402930279415, -3.850878800349884390),
            /*[ 4]*/ (-4.163545402930279415, +3.850878800349884390)
        );

        var z_zeros = ToZArray(high_pass_zeros, dt);
        var z_poles = ToZArray(high_pass_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/  1,
            /*[ 1]*/ (0.928328869671302237, -0.371760016320747899),
            /*[ 2]*/ (0.928328869671302237, +0.371760016320747899),
            /*[ 3]*/ (0.972004020548271175, -0.234964218633382821),
            /*[ 4]*/ (0.972004020548271175, +0.234964218633382821)
        );


        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/  0.590683358317104057,
            /*[ 1]*/ (0.710240523847416094, -0.493566125901245822),
            /*[ 2]*/ (0.710240523847416094, +0.493566125901245822),
            /*[ 3]*/ (0.614384193523889177, -0.257280037461190936),
            /*[ 4]*/ (0.614384193523889177, +0.257280037461190936)
        );

        var k_zeros = z_zeros.Multiply(z => 1 + z).Re;
        var k_poles = z_poles.Multiply(z => 1 + z).Re;

        var g_norm = N.IsEven()
            ? Gp * k_poles / k_zeros
            : 1 * k_poles / k_zeros;

        //g_norm.ToDebug();

        //var zz0 = Complex.Exp(-Consts.pi2 * 0 / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fs / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fp / fd);
        //var zz0 = Complex.Exp(-Consts.pi2 * fd / 2 / fd);
        //var P0 = (1 - zz0).Power(N);
        //var Pp = z_poles.Multiply(z => 1 - z * zz0);
        //var k0 = g_norm * P0 / Pp;
        //k0.Abs.ToDebug();

        //var B = GetCoefficientsInverted(z_zeros).ToArray(z => z.Re * g_norm);

        var B = GetCoefficientsInverted(z_zeros).ToArray(z => z.Re * g_norm);
        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        B.AssertEquals(
            /*[ 0]*/ +0.4427605613526005,
            /*[ 1]*/ -2.1255454758134564,
            /*[ 2]*/ +4.1663872231272000,
            /*[ 3]*/ -4.1663872231272010,
            /*[ 4]*/ +2.1255454758134573,
            /*[ 5]*/ -0.4427605613526008
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ +1,
            /*[ 1]*/ -3.2399327930597144,
            /*[ 2]*/ +4.5020198362216200,
            /*[ 3]*/ -3.2843180660459663,
            /*[ 4]*/ +1.2470797212785330,
            /*[ 5]*/ -0.19603610398067925
        );

        //Debug.WriteLine("---------------------------------------------");

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        filter.B.AssertEquals(Accuracy.Eps(1e-14), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-14), A);

        //var h_0 = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, 0, dt);
        //var h_p = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fs, dt);
        //var h_s = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fp, dt);
        //var h_d = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fd / 2, dt);

        //filter.GetTransmissionCoefficient(0).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fs).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fp).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fd / 2).Abs.ToDebug();

        //var h_0 = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, 0, dt);
        //var h_p = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fs, dt);
        //var h_s = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fp, dt);
        //var h_d = DoubleArrayDSPExtensions.GetTransmissionCoefficient(A, B, fd / 2, dt);

        //h_0.Abs.ToDebug();
        //h_p.Abs.ToDebug();
        //h_s.Abs.ToDebug();
        //h_d.Abs.ToDebug();
    }

    [TestMethod, Ignore]
    public void TransmissionCoefficient_Even()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs);

        var transmission__0 = filter.GetTransmissionCoefficient(0);
        var transmission_fs = filter.GetTransmissionCoefficient(fs);
        var transmission_fp = filter.GetTransmissionCoefficient(fp);
        var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs, 4e-13);
        transmission_fs.Power.In_dB_byPower().AssertEquals(-Rs, 5e-3);
        transmission_fp.Power.In_dB_byPower().AssertEquals(-Rp, 6e-13);
        transmission_fd05.Power.In_dB_byPower().AssertEquals(-Rp, 1.3e-15);
    }

    [TestMethod, Ignore]
    public void TransmissionCoefficient_Odd()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs);

        var transmission__0 = filter.GetTransmissionCoefficient(0);
        var transmission_fs = filter.GetTransmissionCoefficient(fs);
        var transmission_fp = filter.GetTransmissionCoefficient(fp);
        var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fs.Power.In_dB_byPower().AssertEquals(-Rs, 6.81e-1);
        transmission_fp.Power.In_dB_byPower().AssertEquals(-Rp, 8.9e-12);
        transmission_fd05.Power.In_dB_byPower().AssertEquals(0, 2e-15);
    }

    [TestMethod, Ignore]
    public void SignalProcessing_Even()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs);

        static EnumerableSignal GetSinSignal(double f0) => MathEnumerableSignal.Sin(dt, f0, (int)(100 * fd / (fs * 0.1)));

        var x_low = GetSinSignal(fs * 0.1);
        var x_s99 = GetSinSignal(fs * 0.99);
        var x_s = GetSinSignal(fs);

        var x_sp = GetSinSignal(fs + (fp - fs) * 0.1);
        var x_ps = GetSinSignal(fp - (fp - fs) * 0.1);

        var x_p = GetSinSignal(fp);

        var x_pd = GetSinSignal(fp + (fd / 2 - fp) * 0.1);
        var x_dp = GetSinSignal(fp + (fd / 2 - fp) * 0.9);

        var x_fd05 = GetSinSignal(0.9 * (fd / 2));

        /* ----------------------------------------------- */

        var y_low = filter.ProcessIndividual(x_low);
        var y_s99 = filter.ProcessIndividual(x_s99);
        var y_s = filter.ProcessIndividual(x_s);

        var y_sp = filter.ProcessIndividual(x_sp);
        var y_ps = filter.ProcessIndividual(x_ps);

        var y_p = filter.ProcessIndividual(x_p);

        var y_pd = filter.ProcessIndividual(x_pd);
        var y_dp = filter.ProcessIndividual(x_dp);

        var y_fd05 = filter.ProcessIndividual(x_fd05);

        /* ----------------------------------------------- */

        var k_low = y_low.Power / x_low.Power;
        var k_s99 = y_s99.Power / x_s99.Power;
        var k_s = y_s.Power / x_s.Power;

        var k_sp = y_sp.Power / x_sp.Power;
        var k_ps = y_ps.Power / x_ps.Power;

        var k_p = y_p.Power / x_p.Power;

        var k_pd = y_pd.Power / x_pd.Power;
        var k_dp = y_dp.Power / x_dp.Power;

        var k_fd05 = y_fd05.Power / x_fd05.Power;

        /* ----------------------------------------------- */

        var k_low_db = k_low.In_dB_byPower();
        var k_s99_db = k_s99.In_dB_byPower();

        var k_s_db = k_s.In_dB_byPower();

        var k_sp_db = k_sp.In_dB_byPower();
        var k_ps_db = k_ps.In_dB_byPower();

        var k_p_db = k_p.In_dB_byPower();

        var k_pd_db = k_pd.In_dB_byPower();
        var k_dp_db = k_dp.In_dB_byPower();

        var k_fd05_db = k_fd05.In_dB_byPower();

        /* --------------------------------------*/ //   -Rs  -Rp
        // ---+----+->
        k_low_db.AssertLessThan(-Rs);               // \  |    |
        k_s99_db.AssertLessThan(-Rs, 2.7e-1);       //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 3e-1);   //   \+    | fs
        //    \    |
        k_sp_db.AssertGreaterThan(-Rs, 0.5);        //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertGreaterThan(-Rs);             //       \ |
        k_ps_db.AssertLessThan(-Rp);                //        \+ fp
        //         \0дБ
        k_p_db.AssertEquals(-Rp, 5.2e-3);           //          |
        k_pd_db.AssertGreaterOrEqualsThan(-Rp);     //          |
        k_dp_db.AssertGreaterOrEqualsThan(-Rp);     //          |
        k_fd05_db.AssertGreaterOrEqualsThan(-Rp);   //          |fd

        /* ----------------------------------------------- */

        var x =
            x_low +
            x_s99 +
            x_s +
            x_sp +
            x_ps +
            x_p +
            x_pd +
            x_dp +
            x_fd05;

        var y = filter.ProcessIndividual(x);

        var X = x.GetSpectrum();
        var Y = y.GetSpectrum();

        var H = Y / X;

        var h_low = H.GetValue(fs * 0.1).Power.In_dB_byPower();
        var h_s99 = H.GetValue(fs * 0.99).Power.In_dB_byPower();
        var h_s = H.GetValue(fs).Power.In_dB_byPower();

        var h_sp = H.GetValue(fs + (fp - fs) * 0.1).Power.In_dB_byPower();
        var h_ps = H.GetValue(fp - (fp - fs) * 0.1).Power.In_dB_byPower();

        var h_p = H.GetValue(fp).Power.In_dB_byPower();
        var h_pd = H.GetValue(fp + (fd / 2 - fp) * 0.1).Power.In_dB_byPower();
        var h_dp = H.GetValue(fp + (fd / 2 - fp) * 0.9).Power.In_dB_byPower();

        var h_fd05 = H.GetValue(0.9 * (fd / 2)).Power.In_dB_byPower();

        /* --------------------------------------*/ //   -Rs  -Rp
        // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
        //    \    |
        h_sp.AssertGreaterThan(-Rs, 1);             //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertGreaterThan(-Rs);                //       \ |
        h_ps.AssertLessThan(-Rp);                   //        \+ fp
        //         \0дБ
        h_p.AssertEquals(-Rp, 8.1e-5);              //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }

    [TestMethod, Ignore]
    public void SignalProcessing_Odd()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs);

        static EnumerableSignal GetSinSignal(double f0) => MathEnumerableSignal.Sin(dt, f0, (int)(100 * fd / (fs * 0.1)));

        var x_low = GetSinSignal(fs * 0.1);
        var x_s99 = GetSinSignal(fs * 0.99);
        var x_s = GetSinSignal(fs);

        var x_sp = GetSinSignal(fs + (fp - fs) * 0.1);
        var x_ps = GetSinSignal(fp - (fp - fs) * 0.1);

        var x_p = GetSinSignal(fp);

        var x_pd = GetSinSignal(fp + (fd / 2 - fp) * 0.1);
        var x_dp = GetSinSignal(fp + (fd / 2 - fp) * 0.9);

        var x_fd05 = GetSinSignal(0.9 * (fd / 2));

        /* ----------------------------------------------- */

        var y_low = filter.ProcessIndividual(x_low);
        var y_s99 = filter.ProcessIndividual(x_s99);
        var y_s = filter.ProcessIndividual(x_s);

        var y_sp = filter.ProcessIndividual(x_sp);
        var y_ps = filter.ProcessIndividual(x_ps);

        var y_p = filter.ProcessIndividual(x_p);

        var y_pd = filter.ProcessIndividual(x_pd);
        var y_dp = filter.ProcessIndividual(x_dp);

        var y_fd05 = filter.ProcessIndividual(x_fd05);

        /* ----------------------------------------------- */

        var k_low = y_low.Power / x_low.Power;
        var k_s99 = y_s99.Power / x_s99.Power;
        var k_s = y_s.Power / x_s.Power;

        var k_sp = y_sp.Power / x_sp.Power;
        var k_ps = y_ps.Power / x_ps.Power;

        var k_p = y_p.Power / x_p.Power;

        var k_pd = y_pd.Power / x_pd.Power;
        var k_dp = y_dp.Power / x_dp.Power;

        var k_fd05 = y_fd05.Power / x_fd05.Power;

        /* ----------------------------------------------- */

        var k_low_db = k_low.In_dB_byPower();
        var k_s99_db = k_s99.In_dB_byPower();

        var k_s_db = k_s.In_dB_byPower();

        var k_sp_db = k_sp.In_dB_byPower();
        var k_ps_db = k_ps.In_dB_byPower();

        var k_p_db = k_p.In_dB_byPower();

        var k_pd_db = k_pd.In_dB_byPower();
        var k_dp_db = k_dp.In_dB_byPower();

        var k_fd05_db = k_fd05.In_dB_byPower();

        /* --------------------------------------*/ //   -Rs  -Rp
        // ---+----+->
        k_low_db.AssertLessThan(-Rs);               // \  |    |
        k_s99_db.AssertLessThan(-Rs, 1.24);         //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 1.05);   //   \+    | fs
        //    \    |
        k_sp_db.AssertGreaterThan(-Rs, 0.5);        //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertGreaterThan(-Rs);             //       \ |
        k_ps_db.AssertLessThan(-Rp);                //        \+ fp
        //         \0дБ
        k_p_db.AssertEquals(-Rp, 1.6e-2);           //          |
        k_pd_db.AssertGreaterOrEqualsThan(-Rp);     //          |
        k_dp_db.AssertGreaterOrEqualsThan(-Rp);     //          |
        k_fd05_db.AssertGreaterOrEqualsThan(-Rp);   //          |fd

        /* ----------------------------------------------- */

        var x =
            x_low +
            x_s99 +
            x_s +
            x_sp +
            x_ps +
            x_p +
            x_pd +
            x_dp +
            x_fd05;

        var y = filter.ProcessIndividual(x);

        var X = x.GetSpectrum();
        var Y = y.GetSpectrum();

        var H = Y / X;

        var h_low = H.GetValue(fs * 0.1).Power.In_dB_byPower();
        var h_s99 = H.GetValue(fs * 0.99).Power.In_dB_byPower();
        var h_s = H.GetValue(fs).Power.In_dB_byPower();

        var h_sp = H.GetValue(fs + (fp - fs) * 0.1).Power.In_dB_byPower();
        var h_ps = H.GetValue(fp - (fp - fs) * 0.1).Power.In_dB_byPower();

        var h_p = H.GetValue(fp).Power.In_dB_byPower();
        var h_pd = H.GetValue(fp + (fd / 2 - fp) * 0.1).Power.In_dB_byPower();
        var h_dp = H.GetValue(fp + (fd / 2 - fp) * 0.9).Power.In_dB_byPower();

        var h_fd05 = H.GetValue(0.9 * (fd / 2)).Power.In_dB_byPower();

        /* --------------------------------------*/ //   -Rs  -Rp
        // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
        //    \    |
        h_sp.AssertLessOrEqualsThan(-Rs);           //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertGreaterThan(-Rs);                //       \ |
        h_ps.AssertLessThan(-Rp);                   //        \+ fp
        //         \0дБ
        h_p.AssertEquals(-Rp, 8.1e-5);              //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }
}