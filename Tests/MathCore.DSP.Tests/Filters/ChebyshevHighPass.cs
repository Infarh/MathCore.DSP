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
    public void TypeI_Odd_Creation()
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
    public void TypeII_Even_Creation()
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

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(0.662286545701310625);

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
            (zeros[i - r], zeros[i - r + 1]) = Complex.Conjugate(0, 1 / Cos(th));
        }

        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ +1.0195911582083184,
            /*[ 1]*/ -1.0195911582083184,
            /*[ 2]*/ +1.2026897738700906,
            /*[ 3]*/ -1.2026897738700906,
            /*[ 4]*/ +1.7999524462728311,
            /*[ 5]*/ -1.7999524462728311,
            /*[ 6]*/ +5.125830895483011,
            /*[ 7]*/ -5.125830895483011
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.094555282859327405, 0.819754047366180849),
            /*[ 1]*/ (-0.094555282859327405, -0.819754047366180849),
            /*[ 2]*/ (-0.330093870216435492, 0.851931017330359475),
            /*[ 3]*/ (-0.330093870216435492, -0.851931017330359475),
            /*[ 4]*/ (-0.725907452544711673, 0.836437315774920753),
            /*[ 5]*/ (-0.725907452544711673, -0.836437315774920753),
            /*[ 6]*/ (-1.281657580105953764, 0.439636107810248811),
            /*[ 7]*/ (-1.281657580105953764, -0.439636107810248811)
        );

        var high_pass_zeros = AnalogBasedFilter.TransformToHighPassW(zeros, Wp);
        if (N.IsOdd())
            high_pass_zeros = high_pass_zeros.Prepend(0);
        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        //high_pass_zeros.ToIm().ToDebugEnum();
        high_pass_zeros.ToIm().AssertEquals(
            /*[ 0]*/ -3.976300380338443,
            /*[ 1]*/ 3.976300380338443,
            /*[ 2]*/ -3.370944692684622,
            /*[ 3]*/ 3.370944692684622,
            /*[ 4]*/ -2.252393233259301,
            /*[ 5]*/ 2.252393233259301,
            /*[ 6]*/ -0.7909353220657545,
            /*[ 7]*/ 0.7909353220657545
        );

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/ (-0.562968189283692966, -4.880694528620864503),
            /*[ 1]*/ (-0.562968189283692966, 4.880694528620864503),
            /*[ 2]*/ (-1.603197764351980581, -4.137653030244414332),
            /*[ 3]*/ (-1.603197764351980581, 4.137653030244414332),
            /*[ 4]*/ (-2.399355012821603239, -2.764691365931373213),
            /*[ 5]*/ (-2.399355012821603239, 2.764691365931373213),
            /*[ 6]*/ (-2.830232210796488790, -0.970830503144920365),
            /*[ 7]*/ (-2.830232210796488790, 0.970830503144920365)
        );

        var z_zeros = ToZArray(high_pass_zeros, dt);
        var z_poles = ToZArray(high_pass_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (0.923951189091279490, -0.382510392246812769),
            /*[ 1]*/ (0.923951189091279490, 0.382510392246812769),
            /*[ 2]*/ (0.944753122110103538, -0.327782760777945104),
            /*[ 3]*/ (0.944753122110103538, 0.327782760777945104),
            /*[ 4]*/ (0.974951320727045601, -0.222418349541105842),
            /*[ 5]*/ (0.974951320727045601, 0.222418349541105842),
            /*[ 6]*/ (0.996876990801502605, -0.078970027292264058),
            /*[ 7]*/ (0.996876990801502605, 0.078970027292264058)
        );


        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.841500351533454261, -0.437086738035546984),
            /*[ 1]*/ (0.841500351533454261, 0.437086738035546984),
            /*[ 2]*/ (0.786058658178190983, -0.342083199895519718),
            /*[ 3]*/ (0.786058658178190983, 0.342083199895519718),
            /*[ 4]*/ (0.758969056289433008, -0.217104758600424885),
            /*[ 5]*/ (0.758969056289433008, 0.217104758600424885),
            /*[ 6]*/ (0.748900270001869206, -0.074370059550829551),
            /*[ 7]*/ (0.748900270001869206, 0.074370059550829551)
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
            /*[ 0]*/ 0.4304223642768027,
            /*[ 1]*/ -3.3061022631152137,
            /*[ 2]*/ 11.241902806035364,
            /*[ 3]*/ -22.099005752008757,
            /*[ 4]*/ 27.465567953084914,
            /*[ 5]*/ -22.099005752008754,
            /*[ 6]*/ 11.24190280603536,
            /*[ 7]*/ -3.3061022631152124,
            /*[ 8]*/ 0.4304223642768025
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ 1,
            /*[ 1]*/ -6.270856672005895,
            /*[ 2]*/ 17.559656116060708,
            /*[ 3]*/ -28.623642890349473,
            /*[ 4]*/ 29.6656081479576,
            /*[ 5]*/ -19.99499137690193,
            /*[ 6]*/ 8.551631224252297,
            /*[ 7]*/ -2.120383401305136,
            /*[ 8]*/ 0.2332328168275685
        );

        //Debug.WriteLine("---------------------------------------------");

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        filter.B.AssertEquals(Accuracy.Eps(1e-12), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-12), A);

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
    public void TypeII_Odd_Creation()
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

            (zeros[i - r], zeros[i - r + 1]) = Complex.Conjugate(0, 1 / Cos(th));
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

    [TestMethod]
    public void TypeII_Corrected_Even_Creation()
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

        var beta = arcsh(eps_s) / N;
        beta.AssertEquals(0.662286545701310625);

        var r = N % 2;
        var zeros = new Complex[N - r];
        var poles = new Complex[N];
        if (N.IsOdd())
            poles[0] = -kW / Sinh(beta);
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var th = dth * (i - r + 1);
            var (im, re) = Complex.Trigonometry.Cos(new(th, beta));
            var k = kW / (re * re + im * im);
            (poles[i], poles[i + 1]) = Complex.Conjugate(k * re, k * im);
            (zeros[i - r], zeros[i - r + 1]) = Complex.Conjugate(0, kW / Cos(th));
        }

        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ 1.359454877611091,
            /*[ 1]*/ -1.359454877611091,
            /*[ 2]*/ 1.6035863651601208,
            /*[ 3]*/ -1.6035863651601208,
            /*[ 4]*/ 2.3999365950304417,
            /*[ 5]*/ -2.3999365950304417,
            /*[ 6]*/ 6.834441193977347,
            /*[ 7]*/ -6.834441193977347
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.126073710479103179, 1.093005396488240910),
            /*[ 1]*/ (-0.126073710479103179, -1.093005396488240910),
            /*[ 2]*/ (-0.440125160288580675, 1.135908023107145892),
            /*[ 3]*/ (-0.440125160288580675, -1.135908023107145892),
            /*[ 4]*/ (-0.967876603392948898, 1.115249754366560930),
            /*[ 5]*/ (-0.967876603392948898, -1.115249754366560930),
            /*[ 6]*/ (-1.708876773474605093, 0.586181477080331748),
            /*[ 7]*/ (-1.708876773474605093, -0.586181477080331748)
        );

        var high_pass_zeros = AnalogBasedFilter.TransformToHighPassW(zeros, Wp);
        if (N.IsOdd())
            high_pass_zeros = high_pass_zeros.Prepend(0);
        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        //high_pass_zeros.ToIm().ToDebugEnum();
        high_pass_zeros.ToIm().AssertEquals(
            /*[ 0]*/ -2.9822252852538327,
            /*[ 1]*/ 2.9822252852538327,
            /*[ 2]*/ -2.5282085195134667,
            /*[ 3]*/ 2.5282085195134667,
            /*[ 4]*/ -1.6892949249444758,
            /*[ 5]*/ 1.6892949249444758,
            /*[ 6]*/ -0.593201491549316,
            /*[ 7]*/ 0.593201491549316
        );

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/ (-0.422226141962769697, -3.660520896465648377),
            /*[ 1]*/ (-0.422226141962769697, 3.660520896465648377),
            /*[ 2]*/ (-1.202398323263985436, -3.103239772683310527),
            /*[ 3]*/ (-1.202398323263985436, 3.103239772683310527),
            /*[ 4]*/ (-1.799516259616202429, -2.073518524448529909),
            /*[ 5]*/ (-1.799516259616202429, 2.073518524448529909),
            /*[ 6]*/ (-2.122674158097367147, -0.728122877358690301),
            /*[ 7]*/ (-2.122674158097367147, 0.728122877358690301)
        );

        var z_zeros = ToZArray(high_pass_zeros, dt);
        var z_poles = ToZArray(high_pass_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (0.956498873151655826, -0.291736020504174964),
            /*[ 1]*/ (0.956498873151655826, 0.291736020504174964),
            /*[ 2]*/ (0.968543471027404568, -0.248844418724204763),
            /*[ 3]*/ (0.968543471027404568, 0.248844418724204763),
            /*[ 4]*/ (0.985832488466227619, -0.167732837227792864),
            /*[ 5]*/ (0.985832488466227619, 0.167732837227792864),
            /*[ 6]*/ (0.998242106406312613, -0.059268009899843585),
            /*[ 7]*/ (0.998242106406312613, 0.059268009899843585)
        );


        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.897682170271824176, -0.340144369709866945),
            /*[ 1]*/ (0.897682170271824176, 0.340144369709866945),
            /*[ 2]*/ (0.847012268997008544, -0.270333659730207232),
            /*[ 3]*/ (0.847012268997008544, 0.270333659730207232),
            /*[ 4]*/ (0.818450911993809926, -0.172966757927768472),
            /*[ 5]*/ (0.818450911993809926, 0.172966757927768472),
            /*[ 6]*/ (0.806143148078542482, -0.059445532511234814),
            /*[ 7]*/ (0.806143148078542482, 0.059445532511234814)
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

        var B = GetCoefficientsInverted(z_zeros).ToArray(z => z.Re * g_norm);
        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        B.AssertEquals(
            /*[ 0]*/ 0.5143741219808673,
            /*[ 1]*/ -4.021497186490405,
            /*[ 2]*/ 13.846823485433651,
            /*[ 3]*/ -27.423707858170296,
            /*[ 4]*/ 34.16801515496893,
            /*[ 5]*/ -27.423707858170303,
            /*[ 6]*/ 13.846823485433655,
            /*[ 7]*/ -4.021497186490405,
            /*[ 8]*/ 0.5143741219808672
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ 1,
            /*[ 1]*/ -6.738576998682371,
            /*[ 2]*/ 20.08344800090471,
            /*[ 3]*/ -34.55375211985258,
            /*[ 4]*/ 37.51569787777869,
            /*[ 5]*/ -26.308080602509754,
            /*[ 6]*/ 11.631983318758563,
            /*[ 7]*/ -2.963775420148898,
            /*[ 8]*/ 0.33308741373827444
        );

        Debug.WriteLine("---------------------------------------------");

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

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
    public void TypeII_Corrected_Odd_Creation()
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
            poles[0] = -kW / Sinh(beta);
        for (var (i, dth) = (r, Consts.pi05 / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var th = dth * (i - r + 1);
            var (im, re) = Complex.Trigonometry.Cos(new(th, beta));
            var k = kW / (re * re + im * im);
            (poles[i], poles[i + 1]) = Complex.Conjugate(k * re, k * im);

            (zeros[i - r], zeros[i - r + 1]) = Complex.Conjugate(0, kW / Cos(th));
        }

        //zeros.ToIm().ToDebugEnum();
        zeros.ToIm().AssertEquals(
            /*[ 0]*/ 2.1029244484765344,
            /*[ 1]*/ -2.1029244484765344,
            /*[ 2]*/ 3.4026032334081595,
            /*[ 3]*/ -3.4026032334081595
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/  -1.575540533713846836,
            /*[ 1]*/ (-0.311831190556909887, 1.221740635279751741),
            /*[ 2]*/ (-0.311831190556909887, -1.221740635279751741),
            /*[ 3]*/ (-1.049598957226317575, 0.970778022597644896),
            /*[ 4]*/ (-1.049598957226317575, -0.970778022597644896)
        );

        var high_pass_zeros = AnalogBasedFilter.TransformToHighPassW(zeros, Wp);
        if (N.IsOdd())
            high_pass_zeros = high_pass_zeros.Prepend(0);
        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        //high_pass_zeros.ToIm().ToDebugEnum();
        high_pass_zeros.ToIm().AssertEquals(
            /*[ 0]*/ 0,
            /*[ 1]*/ -1.92788700188945,
            /*[ 2]*/ 1.92788700188945,
            /*[ 3]*/ -1.1914996936368127,
            /*[ 4]*/ 1.1914996936368127
        );

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/  -2.573212572714287827,
            /*[ 1]*/ (-0.795166415107994906, -3.115426392761265539),
            /*[ 2]*/ (-0.795166415107994906, 3.115426392761265539),
            /*[ 3]*/ (-2.081772701465139708, -1.925439400174942195),
            /*[ 4]*/ (-2.081772701465139708, 1.925439400174942195)
        );

        var z_zeros = ToZArray(high_pass_zeros, dt);
        var z_poles = ToZArray(high_pass_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/  1,
            /*[ 1]*/ (0.981587346535519156, -0.191013824424721645),
            /*[ 2]*/ (0.981587346535519156, 0.191013824424721645),
            /*[ 3]*/ (0.992926746641183966, -0.118728580403179032),
            /*[ 4]*/ (0.992926746641183966, 0.118728580403179032)
        );


        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/  0.772011842406100768,
            /*[ 1]*/ (0.881299077197208569, -0.281846688830509118),
            /*[ 2]*/ (0.881299077197208569, 0.281846688830509118),
            /*[ 3]*/ (0.797780032897538738, -0.156758995529333173),
            /*[ 4]*/ (0.797780032897538738, 0.156758995529333173)
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
            /*[ 0]*/ 0.6609827401950336,
            /*[ 1]*/ -3.271222211918332,
            /*[ 2]*/ 6.509097273376015,
            /*[ 3]*/ -6.509097273376013,
            /*[ 4]*/ 3.271222211918331,
            /*[ 5]*/ -0.6609827401950334
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ 1,
            /*[ 1]*/ -4.130170062595595,
            /*[ 2]*/ 6.92202112489826,
            /*[ 3]*/ -5.873536007971762,
            /*[ 4]*/ 2.5199790745804,
            /*[ 5]*/ -0.43689818093273614
        );

        //Debug.WriteLine("---------------------------------------------");

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

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
    public void TypeI_Even_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.000188284860342486, +0.001479956305863380, +0.004749902982439790, +0.006626656748170020, -0.002508700952507990,
            -0.027036972447236900, -0.047658930380903000, -0.026875522222909500, +0.048156906403668100, +0.127518952640903000,
            +0.125471055413688000, +0.007199783352445440, -0.153274544241392000, -0.222985771200800000, -0.133943816158011000,
            +0.045014843347039800, +0.170318274178307000, +0.159405560148488000, +0.058976711036568100, -0.020495460062454300,
            -0.024123613028188900, +0.005381232055974190, -0.002406472562804330, -0.056100501168090500, -0.097131880408097200,
            -0.075123193121255600, -0.007892828803182270, +0.042639905274967700, +0.042734268298674500, +0.018681384049698800,
            +0.015751425787796900, +0.040550010009382200, +0.056801780536469400, +0.034122538892406300, -0.013715963681112100,
            -0.045948160711845900, -0.042076589654004700, -0.020758959260145500, -0.011673209094949300, -0.019473763936439400,
            -0.023539106620881500, -0.008750770874050140, +0.014624168091743000, +0.025583569584433500, +0.019550397773783100,
            +0.012234365326072100, +0.016944025107528000, +0.026391792406308700, +0.022572682676533000, +0.001095596925168210,
            -0.021775997147370400, -0.028115701574587000, -0.018918395689351700, -0.010728988422396600, -0.013716200461413900,
            -0.019425093797297700, -0.013011846718691000, +0.006514122496273820, +0.023827839976293800, +0.025283525949237100,
            +0.014432267686419000, +0.006218445501650210, +0.008174708981836740, +0.012176831137656800, +0.006340224134454350,
            -0.008914763914320470, -0.020771674472653600, -0.019371276856041700, -0.009114767408211420, -0.002390764047915550,
            -0.004377544528578730, -0.007877341326967550, -0.003687979784521760, +0.007285462503163000, +0.015172625325745700,
            +0.013317421814009700, +0.005905700026842920, +0.002033447916076970, +0.004384102869773070, +0.006718852935048130,
            +0.002611809965597700, -0.006110425846211280, -0.011644978161238900, -0.009770327310897160, -0.004510704406457110,
            -0.002520749128567200, -0.004776255912809090, -0.005961753341349320, -0.001782503217825400, +0.005365472713167630,
            +0.009318096150187050, +0.007583322840720870, +0.003897148646025110, +0.003080374478408660, +0.005028654108295740,
            +0.005299092144551790, +0.001146184032625430, -0.004811040624158530, -0.007729208292000710, -0.006268888476856120
        };

        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod, Ignore]
    public void TypeI_Odd_ImpulseResponse()
    {
        double[] expected_h =
        {
            +5.7827877333733e-005, +0.000519645569176767, +0.001940177682683280, +0.003362853166507800, +0.000111463882881050,
            -0.012515732982596400, -0.027381321814588900, -0.020802327978941400, +0.025126310539700000, +0.086255814532957000,
            +0.096299248323145800, +0.009212156002258640, -0.130401330094283000, -0.201284751052029000, -0.117511573544446000,
            +0.070891135096577100, +0.209813497749611000, +0.188041104753759000, +0.045103812949553300, -0.079051174566505000,
            -0.094198557190624300, -0.039874456385998300, -0.013837600936063800, -0.047237945512662400, -0.077791963692760100,
            -0.041494082742196500, +0.041533087480998200, +0.092985182877335500, +0.072955759062066600, +0.021921026797681900,
            +0.001494253099410410, +0.017453921488334300, +0.022740623297132400, -0.012025107661415900, -0.057672325144649000,
            -0.065947200445340500, -0.032261968070538700, +0.001244744617363310, +0.002207099318371930, -0.012912636873671800,
            -0.006105465945587130, +0.028808419071681800, +0.057615373634768100, +0.049943767338122200, +0.016057287710447100,
            -0.009450032236347850, -0.011135452252837700, -0.005909641865932880, -0.013696172496154800, -0.029273929112872200,
            -0.032111887189457200, -0.016603545873566100, +0.000416640158841018, +0.003413648982314690, -0.001569922295619740,
            +0.002300443889987560, +0.017514208585507400, +0.028238199404819700, +0.022183547489474400, +0.007267287738381420,
            +0.000280675284928631, +0.004021769226522670, +0.004666323108506110, -0.007639556152399440, -0.023883316927296100,
            -0.027429491512297200, -0.015701400458900900, -0.002958012517850590, -0.000983630371410803, -0.004298625913755040,
            +0.000497235250104881, +0.014818029606027200, +0.025171758121561500, +0.020792218661758700, +0.007570543693995500,
            -0.000149653919325934, +0.001750736605729410, +0.003089076325810160, -0.005305308994708100, -0.018108455455086400,
            -0.022047333565259600, -0.013267551943819000, -0.002018363284070140, +0.001069639556348870, -0.001860209176511530,
            -0.000463614071435545, +0.008704543543861430, +0.017143085286423300, +0.015638408637062000, +0.006296408560066590,
            -0.000535925742056348, +0.000171657836788521, +0.002440508360430670, -0.001418894179803900, -0.009776930284869270,
            -0.013854882455851700, -0.009509493715276500, -0.002402759588650490, -0.000126769418915459, -0.002455878063427290
        };

        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.I);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(5e-11);
    }

    [TestMethod, Ignore]
    public void TypeII_Even_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.014978759442673300, +0.024213566041449500, +0.022188255529808500, -0.011118960862351900, -0.052771796615212700,
            -0.080990391362511100, -0.059416835744030200, +0.018009279363915900, +0.107669861136537000, +0.144958782721866000,
            +0.094009385661312700, -0.021122818851714300, -0.128738485213722000, -0.159246506507874000, -0.096252501126483800,
            +0.012003241661359500, +0.091505253379345200, +0.099587968417279500, +0.052340878657204700, +0.001509352218966960,
            -0.014502171505733300, -0.000066915739404673, +0.011245496870296200, -0.004866090952279130, -0.039066233202330200,
            -0.059562900127856200, -0.044524429264181300, -0.002872839695333230, +0.035092115121871100, +0.045591979860382900,
            +0.030194823189362100, +0.009613926996574650, +0.001633331520038740, +0.005755542436792900, +0.007795207654586580,
            -0.002889288328895110, -0.021820368588420900, -0.034075411598412800, -0.029024608294319500, -0.010059478500377000,
            +0.009337592000632010, +0.018035992925089600, +0.015788021548616600, +0.010682339229426800, +0.009611627835045610,
            +0.011764306321616600, +0.010867695971714800, +0.002795543506369170, -0.009563602973777230, -0.018700185992591400,
            -0.019214859862975300, -0.012168906991686600, -0.003324962859813360, +0.002335603712577620, +0.004327316117023350,
            +0.005542960377646000, +0.008194999027986320, +0.011199936585246900, +0.011396447010625300, +0.006978848930884780,
            -0.000441605742923806, -0.007050968121008230, -0.009929744349876560, -0.009039938905440080, -0.006499090561129910,
            -0.004237239095878490, -0.002444494056096580, -9.84605455818806e-05, +0.003309791904375270, +0.006766991835553470,
            +0.008428940290460120, +0.007229624265207330, +0.003808240514366860, -5.45888262092609e-05, -0.002839383880416040,
            -0.004199726127331250, -0.004693247428208400, -0.004826643271776340, -0.004459058883598450, -0.003093618062839660,
            -0.000639021287652974, +0.002200299271199510, +0.004307562306678230, +0.004932093560538300, +0.004143227345328140,
            +0.002608371888095290, +0.000990977849630774, -0.000447910393354501, -0.001750804133498940, -0.002890320597565150,
            -0.003569851122642080, -0.003408708135062860, -0.002307790529783660, -0.000620292178357255, +0.001040233167470430,
            +0.002174606643479950, +0.002631124732550260, +0.002541032735334680, +0.002085216525012540, +0.001342755609884550
        };

        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod, Ignore]
    public void TypeII_Odd_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.018631883268967100, +0.033242378741733600, +0.025793939719573800, -0.014547489430958800, -0.080881877450630700,
            -0.109077012595877000, -0.063683051969402300, +0.049501110342035800, +0.162026375489273000, +0.182326480269357000,
            +0.078879220643182900, -0.081909952641598600, -0.186256605004734000, -0.166140869133506000, -0.051553968352943200,
            +0.061993882347085800, +0.098313173846301800, +0.059870858088622200, +0.010291414015800200, +0.003314858664158320,
            +0.031686960128460900, +0.044649805263802800, +0.009512302725429230, -0.051064358540042100, -0.083945787411508100,
            -0.061860810584356800, -0.009315765694124060, +0.026696622639163200, +0.025452330103701100, +0.007823890519983960,
            +0.006137645801748450, +0.026165815297591900, +0.043118328277092200, +0.032317623393412800, -0.002385758583597520,
            -0.033245557694446500, -0.038331887224896800, -0.022210823334034400, -0.006786818929099280, -0.005678486468247660,
            -0.011477541197403700, -0.007865733819720940, +0.009689623369064270, +0.028384580968652200, +0.032341152421605100,
            +0.019343345055753400, +0.001831206077536050, -0.007215764802567310, -0.006821477678592540, -0.005993357715613420,
            -0.011344485500707200, -0.019183953313584500, -0.020272046438802700, -0.010704088953410000, +0.003610915862057740,
            +0.013182598752932500, +0.014214132873217800, +0.010857039389329100, +0.009051521629415150, +0.009689654101221560,
            +0.008561244882763530, +0.002261697895311400, -0.007175215814779680, -0.013993186999123500, -0.014421020126974700,
            -0.009936150853576890, -0.004942614199611080, -0.002008600741013500, +1.25620621392513e-05, +0.003440708469094610,
            +0.008241285381879060, +0.011605531454022900, +0.010868206169902300, +0.006302061712568230, +0.000782919153464336,
            -0.003032644384844520, -0.004819226572961160, -0.005942479099369400, -0.007216890403585610, -0.007768177760300830,
            -0.006148611010283620, -0.002230698010497550, +0.002320657304416670, +0.005451388077722360, +0.006393031208349630,
            +0.005887510459671300, +0.004969454929373870, +0.003820732405842630, +0.001908532820429470, -0.000949614092876955,
            -0.003939000233013840, -0.005769546869653470, -0.005761086233614030, -0.004322309109938610, -0.002403241617607810,
            -0.000614444663324253, +0.001050724975163450, +0.002736020379905210, +0.004171901667848610, +0.004729646654726940
        };

        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.II);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod, Ignore]
    public void TypeII_Corrected_Even_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.027988053972296900, +0.062810840221972700, +0.057212037751087700, -0.032413524924819800, -0.155449712245133000,
            -0.201340089046208000, -0.081673207212271200, +0.122360788589317000, +0.239648442095990000, +0.178481468937220000,
            +0.016820632799034200, -0.086866172714557300, -0.068732968535150000, -0.010618063442382300, -0.013306747348527100,
            -0.068522705216921600, -0.090065978394327000, -0.039771381388358800, +0.028306148290808300, +0.048769564366496900,
            +0.026073260320450100, +0.012613418832042100, +0.031085012092532300, +0.051014261117942300, +0.039074708787122900,
            +0.004267364515620640, -0.019383164023757200, -0.018791067384440100, -0.012296484961069800, -0.018371916649107000,
            -0.030970880744940000, -0.031907600388687300, -0.017119594401438300, -0.000256477956192067, +0.006811438529769910,
            +0.007294472100851310, +0.010956994535936600, +0.019146544765477500, +0.023977910510257100, +0.020020911111180400,
            +0.011042133278308000, +0.004035087379887450, +0.000554063658137516, -0.003150167303346150, -0.009236415005702550,
            -0.014924187485088200, -0.016402496105707900, -0.013626288742930800, -0.009684048111074260, -0.006572895028375020,
            -0.003433075868223220, +0.001067315443497130, +0.006173371171857250, +0.009730794638366240, +0.010689767152999900,
            +0.009915497403778650, +0.008620546295284320, +0.006886818484718830, +0.004119693529588900, +0.000415584335887277,
            -0.003162436822515260, -0.005616400439738800, -0.006840478789520730, -0.007272927282160740, -0.007078783152482900,
            -0.006010566733276210, -0.003978075030179170, -0.001423568708404700, +0.000988182436638831, +0.002886175854489640,
            +0.004263785598642890, +0.005142645858307640, +0.005370161231978280, +0.004798329326979330, +0.003537462343457760,
            +0.001927166051533770, +0.000288137468662351, -0.001219483644736440, -0.002507330129167190, -0.003432288730843160,
            -0.003828381519410630, -0.003641368697461960, -0.002977827235270750, -0.002015654647824830, -0.000901641226357868,
            +0.000252788740826069, +0.001317387767456570, +0.002136996605183060, +0.002593872953726820, +0.002660705619384410,
            +0.002385877404102190, +0.001842679994235150, +0.001106763792112070, +0.000271443125875640, -0.000543631188453204,
            -0.001218837256794820, -0.001672492737539670, -0.001873069751020860, -0.001823706496663390, -0.001547222208480880
        };

        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
    }

    [TestMethod, Ignore]
    public void TypeII_Corrected_Odd_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.032046244824887300, +0.071627442221126000, +0.060799080936585900, -0.044935648045190500, -0.189632067746830000,
            -0.219497096142217000, -0.066695223546211800, +0.173081196083005000, +0.282179129693172000, +0.171936487965410000,
            -0.022446697385719100, -0.111747602592483000, -0.062110670666474000, -0.000531906728658865, -0.029716653391800800,
            -0.099724695816931000, -0.097090905058943000, -0.010744701295503700, +0.060196803545379100, +0.051199197847827000,
            +0.011518601591762000, +0.016242892467179400, +0.056606091606206000, +0.065228669313933400, +0.021292229174587800,
            -0.024313520835196600, -0.027927517982712900, -0.009972603193427500, -0.013577740972597900, -0.038647834651119500,
            -0.048113085845902900, -0.026260025824378800, +0.001279069010767870, +0.008225346998165540, +0.002504015613617340,
            +0.007483793771032050, +0.024825123672031500, +0.034354734895956400, +0.025361516861845200, +0.010293511911602200,
            +0.004324920942429720, +0.004834465232773820, -0.000537520990923501, -0.013391687773719900, -0.022669214745099000,
            -0.021215176365308400, -0.014664433497086400, -0.011425542952586400, -0.010849699907745400, -0.006421697695639710,
            +0.003062209213705850, +0.011492179457759300, +0.014124150719775600, +0.013162037902239100, +0.012991834722225900,
            +0.013483359754173400, +0.011089646967685900, +0.004897535715734340, -0.001762749497556640, -0.005814039168720590,
            -0.007854105040952720, -0.009977478402225940, -0.012050445384660200, -0.011984075844528000, -0.008999199505855320,
            -0.004792257483439110, -0.001246033112139870, +0.001618899431696250, +0.004700603099045530, +0.007778139663586220,
            +0.009500192209220670, +0.009148701781380410, +0.007442805356026290, +0.005406976783071380, +0.003178188087519690,
            +0.000389169577457936, -0.002728995299207660, -0.005261150013577360, -0.006562661679887330, -0.006798529804880230,
            -0.006410410020473250, -0.005445916864447090, -0.003705962980986550, -0.001342409398865180, +0.001050655635921650,
            +0.002951144100175120, +0.004261218778096470, +0.005074917926704170, +0.005316688006076580, +0.004799212495150650,
            +0.003554144453966710, +0.001912703378550080, +0.000235285032879509, -0.001310070773379400, -0.002657126594610970,
            -0.003654969650061180, -0.004095336801310540, -0.003902707202644900, -0.003202857885700300, -0.002189151222260010
        };

        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevFilter.ChebyshevType.IICorrected);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.AssertLessThan(3e-11);
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