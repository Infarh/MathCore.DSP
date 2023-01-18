using System.Diagnostics;
using System.Text;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

using static System.Math;

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

        var Fp = ToDigitalFrequency(fp, dt);  // Частота пропускания аналогового прототипа
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.I);

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

        var Fp = ToDigitalFrequency(fp, dt);  // Частота пропускания аналогового прототипа
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.I);

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

        var Fp = ToDigitalFrequency(fp, dt);  // Частота пропускания аналогового прототипа
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
            /*[ 6]*/ +5.1258308954830110,
            /*[ 7]*/ -5.1258308954830110
        );

        //poles.ToDebugEnum();
        poles.AssertEquals(
            /*[ 0]*/ (-0.094555282859327405, +0.819754047366180849),
            /*[ 1]*/ (-0.094555282859327405, -0.819754047366180849),
            /*[ 2]*/ (-0.330093870216435492, +0.851931017330359475),
            /*[ 3]*/ (-0.330093870216435492, -0.851931017330359475),
            /*[ 4]*/ (-0.725907452544711673, +0.836437315774920753),
            /*[ 5]*/ (-0.725907452544711673, -0.836437315774920753),
            /*[ 6]*/ (-1.281657580105953764, +0.439636107810248811),
            /*[ 7]*/ (-1.281657580105953764, -0.439636107810248811)
        );

        var high_pass_zeros = AnalogBasedFilter.TransformToHighPassW(zeros, Wp);
        if (N.IsOdd())
            high_pass_zeros = high_pass_zeros.Prepend(0);
        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        //high_pass_zeros.ToIm().ToDebugEnum();
        high_pass_zeros.ToIm().AssertEquals(
            /*[ 0]*/ -3.9763003803384430,
            /*[ 1]*/ +3.9763003803384430,
            /*[ 2]*/ -3.3709446926846220,
            /*[ 3]*/ +3.3709446926846220,
            /*[ 4]*/ -2.2523932332593010,
            /*[ 5]*/ +2.2523932332593010,
            /*[ 6]*/ -0.7909353220657545,
            /*[ 7]*/ +0.7909353220657545
        );

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/ (-0.562968189283692966, -4.880694528620864503),
            /*[ 1]*/ (-0.562968189283692966, +4.880694528620864503),
            /*[ 2]*/ (-1.603197764351980581, -4.137653030244414332),
            /*[ 3]*/ (-1.603197764351980581, +4.137653030244414332),
            /*[ 4]*/ (-2.399355012821603239, -2.764691365931373213),
            /*[ 5]*/ (-2.399355012821603239, +2.764691365931373213),
            /*[ 6]*/ (-2.830232210796488790, -0.970830503144920365),
            /*[ 7]*/ (-2.830232210796488790, +0.970830503144920365)
        );

        var z_zeros = ToZArray(high_pass_zeros, dt);
        var z_poles = ToZArray(high_pass_poles, dt);

        //z_zeros.ToDebugEnum();
        z_zeros.AssertEquals(
            /*[ 0]*/ (0.923951189091279490, -0.382510392246812769),
            /*[ 1]*/ (0.923951189091279490, +0.382510392246812769),
            /*[ 2]*/ (0.944753122110103538, -0.327782760777945104),
            /*[ 3]*/ (0.944753122110103538, +0.327782760777945104),
            /*[ 4]*/ (0.974951320727045601, -0.222418349541105842),
            /*[ 5]*/ (0.974951320727045601, +0.222418349541105842),
            /*[ 6]*/ (0.996876990801502605, -0.078970027292264058),
            /*[ 7]*/ (0.996876990801502605, +0.078970027292264058)
        );


        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.841500351533454261, -0.437086738035546984),
            /*[ 1]*/ (0.841500351533454261, +0.437086738035546984),
            /*[ 2]*/ (0.786058658178190983, -0.342083199895519718),
            /*[ 3]*/ (0.786058658178190983, +0.342083199895519718),
            /*[ 4]*/ (0.758969056289433008, -0.217104758600424885),
            /*[ 5]*/ (0.758969056289433008, +0.217104758600424885),
            /*[ 6]*/ (0.748900270001869206, -0.074370059550829551),
            /*[ 7]*/ (0.748900270001869206, +0.074370059550829551)
        );

        var k_zeros = z_zeros.Multiply(z => 1 + z).Re;
        var k_poles = z_poles.Multiply(z => 1 + z).Re;

        var g_norm = k_poles / k_zeros;

        g_norm.ToDebug();
        g_norm.AssertEquals(0.48294183586285494, 1.76e-16);

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
            /*[ 0]*/ 0.4829418358628548,
            /*[ 1]*/ -3.709507751024755,
            /*[ 2]*/ 12.613622409840705,
            /*[ 3]*/ -24.795492275479063,
            /*[ 4]*/ 30.816874101245876,
            /*[ 5]*/ -24.79549227547906,
            /*[ 6]*/ 12.613622409840701,
            /*[ 7]*/ -3.7095077510247534,
            /*[ 8]*/ 0.48294183586285455
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.II);

        filter.B.AssertEquals(Accuracy.Eps(1e-12), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-12), A);

        //var h_0 = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, 0, dt);
        //var h_p = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fs, dt);
        //var h_s = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fp, dt);
        //var h_d = g_norm * DoubleArrayDSPExtensions.GetDigitalTransmissionCoefficientFromZPoles(z_zeros, z_poles, fd / 2, dt);

        //var wh = Wp / (ws / wp);
        //var fh = wh / Consts.pi2;

        //filter.GetTransmissionCoefficient(0).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fs).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fp).Abs.ToDebug();
        //filter.GetTransmissionCoefficient(fh).Abs.ToDebug();
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

        var Fp = ToDigitalFrequency(fp, dt);  // Частота пропускания аналогового прототипа
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

        var g_norm = k_poles / k_zeros;

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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.II);

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

        var Fp = ToDigitalFrequency(fp, dt);  // Частота пропускания аналогового прототипа
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

        var g_norm = k_poles / k_zeros;
        g_norm.AssertEquals(0.577137257277902327);
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
            /*[ 0]*/ +00.5771372572779023,
            /*[ 1]*/ -04.5121940571656590,
            /*[ 2]*/ +15.5363914841183920,
            /*[ 3]*/ -30.7699063022528460,
            /*[ 4]*/ +38.3371435507442940,
            /*[ 5]*/ -30.7699063022528530,
            /*[ 6]*/ +15.5363914841183970,
            /*[ 7]*/ -04.5121940571656590,
            /*[ 8]*/ +00.5771372572779022
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ +01,
            /*[ 1]*/ -06.7385769986823710,
            /*[ 2]*/ +20.0834480009047100,
            /*[ 3]*/ -34.5537521198525800,
            /*[ 4]*/ +37.5156978777786900,
            /*[ 5]*/ -26.3080806025097540,
            /*[ 6]*/ +11.6319833187585630,
            /*[ 7]*/ -02.9637754201488980,
            /*[ 8]*/ +00.33308741373827444
        );

        //Debug.WriteLine("---------------------------------------------");

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.IICorrected);

        filter.B.AssertEquals(Accuracy.Eps(1e-13), B);
        filter.A.AssertEquals(Accuracy.Eps(1e-13), A);

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

        var Fp = ToDigitalFrequency(fp, dt);  // Частота пропускания аналогового прототипа
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

        var g_norm = k_poles / k_zeros;

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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.IICorrected);

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
    public void TypeI_Odd_Even_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.I);

        var transmission__0 = filter.FrequencyResponse(0);
        var transmission_fs = filter.FrequencyResponse(fs);
        var transmission_fp = filter.FrequencyResponse(fp);
        var transmission_fd05 = filter.FrequencyResponse(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fs.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fp.Power.In_dB_byPower().AssertEquals(-Rp, 3.6e-10);
        transmission_fd05.Power.In_dB_byPower().AssertEquals(-Rp, 1.8e-15);
    }

    [TestMethod]
    public void TypeI_Odd_Odd_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.I);

        var transmission__0 = filter.FrequencyResponse(0);
        var transmission_fs = filter.FrequencyResponse(fs);
        var transmission_fp = filter.FrequencyResponse(fp);
        var transmission_fd05 = filter.FrequencyResponse(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fs.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fp.Power.In_dB_byPower().AssertEquals(-Rp, 3.6e-10);
        transmission_fd05.Power.In_dB_byPower().AssertEquals(0, 9.7e-16);
    }

    [TestMethod]
    public void TypeII_Odd_Even_TransmissionCoefficient()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 3 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        const double ws = Consts.pi2 * fs;
        const double wp = Consts.pi2 * fp;

        var Fp = ToDigitalFrequency(fp, dt);
        var Wp = Consts.pi2 * Fp;

        var wh = Wp / (ws / wp);
        var fh = wh / Consts.pi2;

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.II);

        var transmission__0 = filter.FrequencyResponse(0);
        var transmission_fs = filter.FrequencyResponse(fs);
        var transmission_fp = filter.FrequencyResponse(fp);
        var transmission_fh = filter.FrequencyResponse(fh);
        var transmission_fd = filter.FrequencyResponse(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fh.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs, 1e-8);
        transmission_fs.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fp.Power.In_dB_byPower().AssertEquals(-Rs, 1.89e-8);
        transmission_fh.Power.In_dB_byPower().AssertGreaterOrEqualsThan(-Rp);
        transmission_fd.Power.In_dB_byPower().AssertEquals(0, 5.8e-15);
    }

    [TestMethod]
    public void TypeII_Odd_Odd_TransmissionCoefficient()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        const double ws = Consts.pi2 * fs;
        const double wp = Consts.pi2 * fp;

        var Fp = ToDigitalFrequency(fp, dt);
        var Wp = Consts.pi2 * Fp;

        var wh = Wp / (ws / wp);
        var fh = wh / Consts.pi2;

        //(fs, fp).ToDebug();

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.II);

        var transmission__0 = filter.FrequencyResponse(0);
        var transmission_fs = filter.FrequencyResponse(fs);
        var transmission_fp = filter.FrequencyResponse(fp);
        var transmission_fh = filter.FrequencyResponse(fh);
        var transmission_fd = filter.FrequencyResponse(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fh.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fs.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fp.Power.In_dB_byPower().AssertEquals(-Rs, 1.89e-8);
        transmission_fh.Power.In_dB_byPower().AssertGreaterOrEqualsThan(-Rp);
        transmission_fd.Power.In_dB_byPower().AssertEquals(0, 5.8e-15);
    }

    [TestMethod]
    public void TypeII_Corrected_Even_Odd_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.IICorrected);

        var transmission__0 = filter.FrequencyResponse(0);
        var transmission_fs = filter.FrequencyResponse(fs);
        var transmission_fp = filter.FrequencyResponse(fp);
        var transmission_fd = filter.FrequencyResponse(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs, 2.1e-7);
        transmission_fs.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fp.Power.In_dB_byPower().AssertGreaterOrEqualsThan(-Rp);
        transmission_fd.Power.In_dB_byPower().AssertEquals(0, 1e-14);
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_Odd_TransmissionCoefficient()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.IICorrected);

        var transmission__0 = filter.FrequencyResponse(0);
        var transmission_fs = filter.FrequencyResponse(fs);
        var transmission_fp = filter.FrequencyResponse(fp);
        var transmission_fd = filter.FrequencyResponse(fd / 2);

        //transmission__0.Abs.In_dB().ToDebug();
        //transmission_fs.Abs.ToDebug();
        //transmission_fp.Abs.ToDebug();
        //transmission_fd05.Abs.ToDebug();

        transmission__0.Power.In_dB_byPower().AssertLessThan(-Rs, 2.1e-7);
        transmission_fs.Power.In_dB_byPower().AssertLessThan(-Rs);
        transmission_fp.Power.In_dB_byPower().AssertGreaterOrEqualsThan(-Rp);
        transmission_fd.Power.In_dB_byPower().AssertEquals(0, 2e-15);
    }

    [TestMethod]
    public void TypeI_Even_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.224811971559488000, -0.632593721592196000, +0.269045081867305000, +0.338963022597473000, +0.071226091837917200,
            -0.115555846809572000, -0.161752313338841000, -0.129072421247452000, -0.061582644281420800, +0.015734064422365000,
            +0.078828799662082100, +0.108044682917050000, +0.097614400425606100, +0.056894523500653300, +0.003886813694557760,
            -0.043060954451063500, -0.071301607409635300, -0.076371614368867400, -0.061267937208092000, -0.033693546754163200,
            -0.002706991102998990, +0.024018761564291600, +0.041689179707460800, +0.048613390962327600, +0.045567373358783800,
            +0.034832000761474000, +0.019329335710893300, +0.002050993428481040, -0.014241080054372500, -0.027175249395771100,
            -0.034924669734368800, -0.036408874136679100, -0.031500415278672100, -0.021136432807733700, -0.007238982604578150,
            +0.007593854481707800, +0.020582105361241800, +0.029329895371685400, +0.032319831429271500, +0.029210176667435400,
            +0.020865346681641200, +0.009123213140285220, -0.003631368755634910, -0.014976557822717500, -0.022928714132305600,
            -0.026286140035436800, -0.024785732216932900, -0.019063817829415500, -0.010453341004051400, -0.000676670854755844,
            +0.008496143122623960, +0.015558342147786200, +0.019496745050476400, +0.019912662464279200, +0.017027035976066600,
            +0.011583050116253800, +0.004674056929378450, -0.002467057249889220, -0.008679188778345610, -0.013051360846780900,
            -0.015043340533055700, -0.014535565245356300, -0.011808001704905800, -0.007458871400068930, -0.002282484376124620,
            +0.002870321567070050, +0.007225705813897420, +0.010194719656387300, +0.011442430480772200, +0.010912741156831400,
            +0.008810686210277810, +0.005549872429671890, +0.001676518082883410, -0.002216935334958380, -0.005576160484510380,
            -0.007956748851913920, -0.009075893762799030, -0.008839600266521940, -0.007344237343137010, -0.004854085614859820,
            -0.001758909207861150, +0.001482570124068470, +0.004406144045461270, +0.006604748378670040, +0.007782561124732490,
            +0.007793767731849430, +0.006660901581628520, +0.004570105079575940, +0.001843607092202990, -0.001107168513063940,
            -0.003841595807511310, -0.005954546131888210, -0.007136393783901750, -0.007219298417376810, -0.006202522762030570,
            -0.004252553082627180, -0.001677596755372440, +0.001119951431237760, +0.003705912951975380, +0.005681928809655460
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.I);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.ToDebug().AssertLessThan(1e-20);
    }

    [TestMethod]
    public void TypeI_Odd_ImpulseResponse()
    {
        double[] expected_h =
        {
            +0.414203053702161000, -0.692809274016991000, -0.065641613135659300, +0.135612285710467000, +0.189913794868845000,
            +0.161125935352069000, +0.087552015077432400, +0.003804032739121860, -0.063456225756700300, -0.099767750876810300,
            -0.103025252550870000, -0.080320752719879400, -0.043337141300307800, -0.004009133337524600, +0.028475539033938300,
            +0.048996197515273200, +0.056378157578725000, +0.052374845673928100, +0.040341236165606500, +0.024048589214155100,
            +0.006876035136656080, -0.008585325027516820, -0.020611252160978200, -0.028249860373842300, -0.031179296356041300,
            -0.029589796482855700, -0.024105332241039600, -0.015721218102915400, -0.005721752701336400, +0.004447186311132420,
            +0.013352102772325700, +0.019769675016485300, +0.022870158170766700, +0.022345885053875400, +0.018454704731131900,
            +0.011966859278468400, +0.004024914521522090, -0.004055823797074380, -0.011008249845472000, -0.015815589520960500,
            -0.017860706282032400, -0.017000527859017400, -0.013558373626919800, -0.008241534186469450, -0.002003258359893080,
            +0.004124305454975220, +0.009197346776399340, +0.012501428516078100, +0.013646325450091300, +0.012603520887082100,
            +0.009684951861291060, +0.005470717190053360, +0.000700667540866244, -0.003850865547607030, -0.007496779356584160,
            -0.009736949877419990, -0.010322301170226200, -0.009275850321513320, -0.006870474769631860, -0.003569787018181010,
            +5.62933889858737e-05, +0.003426570842243710, +0.006038877272249220, +0.007542071394815400, +0.007779885960436650,
            +0.006801771212912290, +0.004841294383463460, +0.002267618835526590, -0.000480633676435929, -0.002967744628583860,
            -0.004827900596202400, -0.005816862770341290, -0.005841373269416320, -0.004963273249710580, -0.003379415677085170,
            -0.001382117264834380, +0.000692425599481376, +0.002518980115467440, +0.003832697180796690, +0.004465853766076140,
            +0.004367112740303660, +0.003601519442648250, +0.002332547038553500, +0.000790182318406353, -0.000769075054517568,
            -0.002103579518850910, -0.003022799400860810, -0.003413224612300030, -0.003250631134697020, -0.002597739482737430,
            -0.001588622616624900, -0.000403170161016594, +0.000763775686575668, +0.001733489556908560, +0.002369831678042250,
            +0.002597346775040560, +0.002408771339986930, +0.001861512328300180, +0.001064391063165980, +0.000157359511452362
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.I);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.ToDebug().AssertLessThan(1e-18);
    }

    private static void ReadImpulseResponse(string FileName)
    {
        var file_path = Path.Combine("Filters", "MathCad", FileName);
        var file = new FileInfo(file_path);
        file.ThrowIfNotFound();

        var lines = file.ReadLines()
           .Select(s => s!.StartsWith('-') ? s : $"+{s}");

        var max_str_len = lines.Max(s => s.Length);

        const int columns_count = 4;
        var values = lines
           .Select(s => s.Length == max_str_len ? s : $"{s}{new string('0', max_str_len - s.Length)}")
           .SelectGroup(columns_count)
           .Select(g => g.JoinStrings(", "))
           .JoinStrings(",\r\n            ");

        var h = new StringBuilder("        double[] expected_h =\r\n")
           .AppendLine("            {")
           .AppendLine($"            {values}")
           .AppendLine("        };");
        h.ToDebug("");
    }

    [TestMethod]
    public void TypeII_Even_ImpulseResponse()
    {
        //ReadImpulseResponse("ChebyshevII.Even.HighPass.Impulse.csv");

        double[] expected_h =
        {
            +0.482941835853120000, -0.681048717443790000, -0.137429045565316000, +0.125245800043701000,
            +0.194620802937324000, +0.152096008083167000, +0.064323856664702900, -0.021425378968517800,
            -0.078112176487666200, -0.096619424582816400, -0.081287214752554700, -0.044449474554657500,
            -0.001079708631932400, +0.035571674230683000, +0.056802470799113600, +0.059575155374870400,
            +0.046096042377029900, +0.022282233372760300, -0.004326321033392610, -0.026647900026603200,
            -0.039676545371301600, -0.041357295367898700, -0.032645529480487800, -0.016864450070844300,
            +0.001390455056124520, +0.017489052825038300, +0.027848615104991700, +0.030637911369313700,
            +0.025994099399247900, +0.015765043981264300, +0.002893125858468670, -0.009378298389328620,
            -0.018277153351721100, -0.022083136153851800, -0.020408181942377700, -0.014147417792351300,
            -0.005147386463829060, +0.004304256579633330, +0.012035204288903400, +0.016462745381322600,
            +0.016889186955833600, +0.013571417288639900, +0.007572219906881950, +0.000448188035939360,
            -0.006140618044863840, -0.010805618594689200, -0.012708477947002600, -0.011692110670078800,
            -0.008249867831198270, -0.003355617642908590, +0.001794461991387730, +0.006063370359745280,
            +0.008614440391736430, +0.009063661601311010, +0.007518811755807100, +0.004508064067912650,
            +0.000824479507822859, -0.002671571118495740, -0.005245062914779710, -0.006432798650323930,
            -0.006116576701793710, -0.004514490002352240, -0.002100382075418590, +0.000523961285422105,
            +0.002771530581792130, +0.004195373090816010, +0.004571092695813570, +0.003923014048901150,
            +0.002493995013649800, +0.000671056756794874, -0.001112667434688570, -0.002476100792500210,
            -0.003167284850477910, -0.003104803257490180, -0.002378198674456820, -0.001211415826266940,
            +9.91111710360388e-05, +0.001255798064924110, +0.002024312767673810, +0.002277823258589060,
            +0.002013559814582430, +0.001340901767827100, +0.000446432330683897, -0.000454165919756734,
            -0.001165653146000060, -0.001553360738780460, -0.001566183082818170, -0.001239191563229830,
            -0.000677360618082777, -2.58276484150249e-05, +0.000565527598104008, +0.000974955806258183,
            +0.001132315034460690, +0.001029025148209250, +0.000713711980832986, +0.000275925975033635,
            -0.000177342771820655, -0.000546550272082854, -0.000760364838253820, -0.000788240484661646
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.II);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.ToDebug().AssertLessThan(1e-20);
    }

    [TestMethod]
    public void TypeII_Odd_ImpulseResponse()
    {
        //ReadImpulseResponse("ChebyshevII.Odd.HighPass.Impulse.csv");

        double[] expected_h =
        {
            +0.442760561352705000, -0.691031013613667000, -0.065823648807067700, +0.185550419860105000, //  1
            +0.201332378851880000, +0.106573164710804000, -0.005085565979994560, -0.079331801927750000, //  2
            -0.098817162964092100, -0.073147054494754100, -0.025431288829409900, +0.020302956985340200, //  3
            +0.047715173023512200, +0.051513674075933000, +0.036142935914436200, +0.011592048846027800, //  4
            -0.011495976903317500, -0.025616690303070500, -0.028143846851301300, -0.020984656724823500, //  5
            -0.008809220771015720, +0.003191142300146940, +0.011153910147461100, +0.013491380246851600, //  6
            +0.010848817093900400, +0.005337431214174860, -0.000522923508235848, -0.004730749874999680, //  7
            -0.006327820062628130, -0.005450671544804500, -0.003010668627598690, -0.000200796101462262, //  8
            +0.001965683251320090, +0.002941617106652970, +0.002707633840717510, +0.001645471119961030, //  9
            +0.000311873781631215, -0.000787848945513283, -0.001352361157815850, -0.001331594663795880, // 10
            -0.000879826706238848, -0.000253649025504895, +0.000297866415261269, +0.000612862454357990, // 11
            +0.000647741990175660, +0.000461652394284696, +0.000171223140709275, -0.000102121814741882, // 12
            -0.000273147226311898, -0.000311609488884479, -0.000238308461513276, -0.000105413707598904, // 13
            +2.85282230595596e-05, +0.000119378499080283, +0.000148234483175034, +0.000121263385359278, // 14
            +6.13657353787184e-05, -3.54206042626946e-06, -5.09359904863435e-05, -6.97039645423488e-05, // 15
            -6.0910536577185e-050, -3.43802762265477e-05, -3.2724613970119e-060, +2.10696623975816e-05, // 16
            +3.23770152806496e-05, +3.02297924111329e-05, +1.8721248436943e-050, +3.97980856798954e-06, // 17
            -8.35117084228995e-06, -1.4839772091063e-050, -1.48325288269083e-05, -9.96844745463066e-06, // 18
            -3.06453040818622e-06, +3.10378106933984e-06, +6.70128150196362e-06, +7.19726267629708e-06, // 19
            +5.21068748058263e-06, +2.0177767787338e-060, -1.03163653061742e-06, -2.97481629466015e-06, // 20
            -3.45395536897492e-06, -2.6809694252563e-060, -1.22453797723041e-06, +2.66063836312821e-07, // 21
            +1.29395319600496e-06, +1.63900717262793e-06, +1.36023807430245e-06, +7.06134526577264e-07, // 22
            -1.44741852154404e-08, -5.48783618195211e-07, -7.687093275448e-0700, -6.81419927089322e-07, // 23
            -3.9291175350389e-070, -4.83881213857461e-08, +2.25182881539987e-07, +3.56064709640052e-07, // 24
            +3.37335314513599e-07, +2.12824525296398e-07, +4.99698209161341e-08, -8.82220162103009e-08, // 25
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.II);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.ToDebug().AssertLessThan(3e-18);
    }

    [TestMethod]
    public void TypeII_Corrected_Even_ImpulseResponse()
    {
        //ReadImpulseResponse("ChebyshevII.Corrected.Even.HighPassImpulse.csv");

        double[] expected_h =
        {
            +0.577137257237147000, -0.623110210255885000, -0.253390741919103000, -0.020940094629099100,
            +0.102494174254646000, +0.145494874375050000, +0.134870572707153000, +0.094032220251121900,
            +0.041915131400757900, -0.007543532406946140, -0.045410778230985500, -0.067326647733892900,
            -0.072743960888151600, -0.063962429805785000, -0.045082737758049500, -0.020996067590109100,
            +0.003495516948070810, +0.024359452847013800, +0.038775945756840000, +0.045336759614385300,
            +0.044019528534572100, +0.035979292514061700, +0.023213325892383700, +0.008161644874059880,
            -0.006695922439137060, -0.019193477426268000, -0.027744282271795400, -0.031494413038049000,
            -0.030358751468800100, -0.024948878900685400, -0.016415041988196000, -0.006232757388672710,
            +0.004031740183112760, +0.012944052279841900, +0.019384128286531100, +0.022668800837617600,
            +0.022604209490705800, +0.019468067773740600, +0.013930819256153600, +0.006931774201129430,
            -0.000469371067836504, -0.007243939166612000, -0.012528517315827100, -0.015725996760281400,
            -0.016563088905985400, -0.015100542685620200, -0.011698735851964200, -0.006946854626162490,
            -0.001567884301962040, +0.003686270312292180, +0.008135134718298090, +0.011252858083437100,
            +0.012724867998949800, +0.012472783296839000, +0.010646883651815400, +0.007589658974420210,
            +0.003777392099743120, -0.000251007613865691, -0.003967798268074760, -0.006923166956083470,
            -0.008796956868512850, -0.009429906587291880, -0.008831981100326430, -0.007168510277138830,
            -0.004727581813105230, -0.001874292157276360, +0.001001325861073140, +0.003535585233252560,
            +0.005435186262722710, +0.006509000066254690, +0.006684483566593440, +0.006007835282073020,
            +0.004629119303806250, +0.002775396210731350, +0.000716183411454840, -0.001273796030019720,
            -0.002949657696828530, -0.004124934765138000, -0.004690580429136460, -0.004622734220101640,
            -0.003979141162431010, -0.002885554832431950, -0.001514617104884790, -6.04637374998533e-05,
            +0.001287407998774830, +0.002367655337023950, +0.003065226663545140, +0.003322249835134980,
            +0.003140863479428640, +0.002578262245313260, +0.001735188177288520, +0.000739831061548968,
            -0.000270473484804181, -0.001167268178040270, -0.001846469293079870, -0.002239576712167080,
            -0.002319520385678260, -0.002100825656071490, -0.001634521293568140, -0.000998844003805030
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.IICorrected);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.ToDebug().AssertLessThan(3e-11);
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_ImpulseResponse()
    {
        //ReadImpulseResponse("ChebyshevII.Corrected.Odd.HighPassImpulse.csv");

        double[] expected_h =
        {
            +0.660982739493786000, -0.541251087174461000, -0.301698251142820000, -0.126304978563997000,
            -0.006797581247875740, +0.065915850230827100, +0.101241497297906000, +0.108422341123950000,
            +0.096113391366116500, +0.072032081697388300, +0.042700407001810000, +0.013288011671961300,
            -0.012442557150599300, -0.032094969191325000, -0.044515868917511000, -0.049607643674429800,
            -0.048098684336954700, -0.041293627956146800, -0.030825053682935000, -0.018425381406094300,
            -0.005733707985777640, +0.005852478612634620, +0.015275817412330500, +0.021867768826259600,
            +0.025351625652738800, +0.025807261042803000, +0.023607212839488900, +0.019334741487808500,
            +0.013694443122191200, +0.007425026248717720, +0.001222159866881990, -0.004322851280571280,
            -0.008765165049771680, -0.011828187898751000, -0.013406029398863500, -0.013549273120488400,
            -0.012437884394850400, -0.010345718071520100, -0.007601220298228930, -0.004548600446267950,
            -0.001513077904166910, +0.001227102174428900, +0.003460220060721870, +0.005051571200049000,
            +0.005945209985313360, +0.006158040072797150, +0.005767905544332510, +0.004897662387165760,
            +0.003697302161182100, +0.002326084401081770, +0.000936352107539958, -0.000339693833833920,
            -0.001399446175144960, -0.002175224411085970, -0.002635580393201140, -0.002783031528042410,
            -0.002648945786916380, -0.002286465389345640, -0.001762412241589110, -0.001149078084372060,
            -0.000516683382058699, +7.28860985623459e-05, +0.000570663182593314, +0.000943316444834068,
            +0.001174015047321720, +0.001261626721308930, +0.001218563586778980, +0.001067671687889620,
            +0.000838592614819978, +0.000564013340296244, +0.000276171001367021, +3.90276732396555e-06,
            -0.000229561772750958, -0.000407966105678969, -0.000522543653421147, -0.000571759774182685,
            -0.000560414530799546, -0.000498282118072610, -0.000398480864749686, -0.000275765021698240,
            -0.000144909350456081, -1.93242412325152e-05, +9.00022875613663e-05, +0.000175183791126110,
            +0.000231726781579536, +0.000258462841783735, +0.000257179512161521, +0.000232027814500384,
            +0.000188793964083477, +0.000134122926603918, +7.47733682036106e-05, +1.69691654566576e-05,
            -3.41057650425988e-05, -7.46426899089079e-05, -0.000102366433366721, -0.000116527534501616,
            -0.000117752835393298, -0.000107788842523960, -8.91772915919568e-05, -6.49030137195547e-05
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.IICorrected);

        var actual_h = filter.GetImpulseResponse(expected_h.Length);

        var error2 = expected_h.Zip(actual_h, (e, a) => (e - a).Pow2()).Average();

        error2.ToDebug().AssertLessThan(3e-11);
    }

    [TestMethod]
    public void TypeI_Even_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.I);

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
        //                                          // ---+----+->
        k_low_db.AssertLessThan(-Rs);               // \  |    |
        k_s99_db.AssertLessThan(-Rs, 2.7e-1);       //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 3e-1);   //   \+    | fs
        //                                          //    \    |
        k_sp_db.AssertGreaterThan(-Rs, 0.5);        //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertGreaterThan(-Rs);             //       \ |
        k_ps_db.AssertLessThan(-Rp);                //        \+ fp
        //                                          //         \0дБ
        k_p_db.AssertEquals(-Rp, 2.19e-2);          //          |
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
        //                                          // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
        //                                          //    \    |
        h_sp.AssertGreaterThan(-Rs, 1);             //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertGreaterThan(-Rs);                //       \ |
        h_ps.AssertLessThan(-Rp);                   //        \+ fp
        //                                          //         \0дБ
        h_p.AssertEquals(-Rp, 8.1e-5);              //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }

    [TestMethod]
    public void TypeI_Odd_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.I);

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
        //                                          // ---+----+->
        k_low_db.AssertLessThan(-Rs);               // \  |    |
        k_s99_db.AssertLessThan(-Rs, 2.7e-1);       //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 3e-1);   //   \+    | fs
        //                                          //    \    |
        k_sp_db.AssertGreaterThan(-Rs, 0.5);        //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertGreaterThan(-Rs);             //       \ |
        k_ps_db.AssertLessThan(-Rp);                //        \+ fp
        //                                          //         \0дБ
        k_p_db.AssertEquals(-Rp, 2.19e-2);          //          |
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
        //                                          // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
        //                                          //    \    |
        h_sp.AssertGreaterThan(-Rs, 1);             //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertGreaterThan(-Rs);                //       \ |
        h_ps.AssertLessThan(-Rp);                   //        \+ fp
        //                                          //         \0дБ
        h_p.AssertEquals(-Rp, 8.1e-5);              //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }

    [TestMethod]
    public void TypeII_Even_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.II);

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
        //                                          // ---+----+->
        k_low_db.AssertLessThan(-Rs);               // \  |    |
        k_s99_db.AssertLessThan(-Rs, 2.7e-1);       //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 3e-1);   //   \+    | fs
        //                                          //    \    |
        k_sp_db.AssertLessThan(-Rs);                //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertLessThan(-Rs);                //       \ |
        k_ps_db.AssertLessThan(-Rp);                //        \+ fp
        //                                          //         \0дБ
        k_p_db.AssertEquals(-Rs, 2.13);             //          |
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
        //                                          // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
        //                                          //    \    |
        h_sp.AssertLessThan(-Rs);                   //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertLessThan(-Rs);                   //       \ |
        h_ps.AssertLessThan(-Rp);                   //        \+ fp
        //                                          //         \0дБ
        h_p.AssertEquals(-Rs, 1.25e-3);             //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }

    [TestMethod]
    public void TypeII_Odd_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.II);

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

        /* --------------------------------------*/  //   -Rs  -Rp
        //                                           // ---+----+->
        k_low_db.AssertLessThan(-Rs);                // \  |    |
        k_s99_db.AssertLessThan(-Rs, 2.7e-1);        //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 3e-1);    //   \+    | fs
        //                                           //    \    |
        k_sp_db.AssertLessThan(-Rs);                 //     \   |
        k_sp_db.AssertLessThan(-Rp);                 //      \  |
        k_ps_db.AssertLessThan(-Rs);                 //       \ |
        k_ps_db.AssertLessThan(-Rp);                 //        \+ fp
        //                                           //         \0дБ
        k_p_db.AssertEquals(-Rs, 2.13);              //          |
        k_pd_db.AssertGreaterOrEqualsThan(-Rp, 0.56);//          |
        k_dp_db.AssertGreaterOrEqualsThan(-Rp);      //          |
        k_fd05_db.AssertGreaterOrEqualsThan(-Rp);    //          |fd

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
        //                                          // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
        //                                          //    \    |
        h_sp.AssertLessThan(-Rs);                   //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertLessThan(-Rs);                   //       \ |
        h_ps.AssertLessThan(-Rp);                   //        \+ fp
        //                                          //         \0дБ
        h_p.AssertEquals(-Rs, 1.25e-3);             //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp, 0.56);  //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }

    [TestMethod]
    public void TypeII_Corrected_Even_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.IICorrected);

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
        //                                          // ---+----+->
        k_low_db.AssertLessThan(-Rs);               // \  |    |
        k_s99_db.AssertLessThan(-Rs, 2.7e-1);       //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 0.75);   //   \+    | fs
        //                                          //    \    |
        k_sp_db.AssertGreaterThan(-Rs, 0.5);        //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertGreaterThan(-Rs);             //       \ |
        k_ps_db.AssertLessThan(-Rp, 0.11);          //        \+ fp
        //                                          //         \0дБ
        k_p_db.AssertEquals(-Rp, 0.51);             //          |
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
        //                                          // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
        //                                          //    \    |
        h_sp.AssertGreaterThan(-Rs, 1);             //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertGreaterThan(-Rs);                //       \ |
        h_ps.AssertLessThan(-Rp, 0.12);             //        \+ fp
        //                                          //         \0дБ
        h_p.AssertEquals(-Rp, 0.52);                //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }

    [TestMethod]
    public void TypeII_Corrected_Odd_SignalProcessing()
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

        var filter = new DSP.Filters.ChebyshevHighPass(dt, fs, fp, Gp, Gs, ChebyshevType.IICorrected);

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
        //                                          // ---+----+->
        k_low_db.AssertLessThan(-Rs);               // \  |    |
        k_s99_db.AssertLessThan(-Rs, 2.7e-1);       //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs, 0.75);   //   \+    | fs
        //                                          //    \    |
        k_sp_db.AssertGreaterThan(-Rs, 0.5);        //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertGreaterThan(-Rs);             //       \ |
        k_ps_db.AssertLessThan(-Rp, 0.43);          //        \+ fp
        //                                          //         \0дБ
        k_p_db.AssertEquals(-Rp, 0.68);             //          |
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
        //                                          // ---+----+->
        h_low.AssertLessThan(-Rs);                  // \  |    |
        h_s99.AssertLessThan(-Rs);                  //  \ |    |
        h_s.AssertLessOrEqualsThan(-Rs);            //   \+    | fs
        //                                          //    \    |
        h_sp.AssertGreaterThan(-Rs, 1);             //     \   |
        h_sp.AssertLessThan(-Rp);                   //      \  |
        h_ps.AssertGreaterThan(-Rs);                //       \ |
        h_ps.AssertLessThan(-Rp, 0.43);             //        \+ fp
        //                                          //         \0дБ
        h_p.AssertEquals(-Rp, 0.69);                //          |
        h_pd.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_dp.AssertGreaterOrEqualsThan(-Rp);        //          |
        h_fd05.AssertGreaterOrEqualsThan(-Rp);      //          |fd
    }
}