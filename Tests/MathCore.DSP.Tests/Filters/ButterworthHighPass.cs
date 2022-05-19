using System;
using System.Diagnostics;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using MathCore.DSP.Tests.Infrastructure;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.VisualStudio.TestTools.UnitTesting.Extensions;

using static System.Math;
using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ButterworthHighPass
{
    [TestMethod]
    public void Creation()
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

        //const double ws = 2 * PI * fs / fd;
        //const double wp = 2 * PI * fp / fd;

        var Fs = DigitalFilter.ToAnalogFrequency(fs, dt);
        var Fp = DigitalFilter.ToAnalogFrequency(fp, dt);

        //(Fp, Fs).ToDebug();

        //var Ws = Consts.pi2 * Fs;
        var Wp = Consts.pi2 * Fp;

        var eps_p = Sqrt(Pow(10, Rp / 10) - 1);
        var eps_s = Sqrt(Pow(10, Rs / 10) - 1);

        eps_p.AssertEquals(0.5088471399095875);
        eps_s.AssertEquals(99.99499987499375);

        var kEps = eps_s / eps_p;
        var kW = Fp / Fs;

        kEps.AssertEquals(196.51284645671973);
        kW.AssertEquals(2.020338844941193);

        var N = (int)Ceiling(Log(kEps) / Log(kW));
        N.AssertEquals(8);

        var alpha = eps_p.Pow(-1d / N);
        alpha.AssertEquals(1.0881194736627366);

        var poles = new Complex[N];
        var r = N % 2;
        if (r != 0) poles[0] = -alpha;
        for (var (i, th0) = (r, Consts.pi05 / N); i < poles.Length; i += 2)
            (poles[i], poles[i + 1]) = Complex.ConjugateAbsExp(alpha, th0 * (i - r + 1 + N));

        //poles.ToDebugEnum();
        poles.AssertEquals(AccuracyComplex.Eps(1e-15),
            /*[ 0]*/ (-0.212281578508883212, +1.067211563088522608),
            /*[ 1]*/ (-0.212281578508883212, -1.067211563088522608),
            /*[ 2]*/ (-0.604526789535973275, +0.904738276905205363),
            /*[ 3]*/ (-0.604526789535973275, -0.904738276905205363),
            /*[ 4]*/ (-0.904738276905205363, +0.604526789535973497),
            /*[ 5]*/ (-0.904738276905205363, -0.604526789535973497),
            /*[ 6]*/ (-1.067211563088522608, +0.212281578508883656),
            /*[ 7]*/ (-1.067211563088522608, -0.212281578508883656)
        );

        var high_pass_poles = AnalogBasedFilter.TransformToHighPassW(poles, Wp);

        //high_pass_poles.ToDebugEnum();
        high_pass_poles.AssertEquals(
            /*[ 0]*/ (-0.726882792937593614, -3.654286571081900803),
            /*[ 1]*/ (-0.726882792937593614, +3.654286571081900803),
            /*[ 2]*/ (-2.069987062796957655, -3.097954566824938194),
            /*[ 3]*/ (-2.069987062796957655, +3.097954566824938194),
            /*[ 4]*/ (-3.097954566824937750, -2.069987062796958099),
            /*[ 5]*/ (-3.097954566824937750, +2.069987062796958099),
            /*[ 6]*/ (-3.654286571081899915, -0.726882792937594946),
            /*[ 7]*/ (-3.654286571081899915, +0.726882792937594946)
        );


        var z_poles = DigitalFilter.ToZArray(high_pass_poles, dt);

        //z_poles.ToDebugEnum();
        z_poles.AssertEquals(
            /*[ 0]*/ (0.871681477170508701, -0.329989827978254768),
            /*[ 1]*/ (0.871681477170508701, +0.329989827978254768),
            /*[ 2]*/ (0.777394997247055186, -0.249492169302408917),
            /*[ 3]*/ (0.777394997247055186, +0.249492169302408917),
            /*[ 4]*/ (0.717957565302678624, -0.153959517251727629),
            /*[ 5]*/ (0.717957565302678624, +0.153959517251727629),
            /*[ 6]*/ (0.689430099459400703, -0.051915227520291297),
            /*[ 7]*/ (0.689430099459400703, +0.051915227520291297)
        );

        var g_norm = z_poles.Multiply(z => (1 + z) / 2).Abs;
        g_norm.AssertEquals(0.3863217046553797);

        var B = new double[N + 1];
        for (var i = 0; i < B.Length; i++)
            B[i] = BinomialCoefficient(N, i) * (i % 2 == 0 ? g_norm : -g_norm);
        var A = GetCoefficientsInverted(z_poles).ToRe();

        //B.ToDebugEnum();
        B.AssertEquals(
            /*[ 0]*/ +00.3863217046553797,
            /*[ 1]*/ -03.0905736372430375,
            /*[ 2]*/ +10.8170077303506320,
            /*[ 3]*/ -21.6340154607012640,
            /*[ 4]*/ +27.0425193258765800,
            /*[ 5]*/ -21.6340154607012640,
            /*[ 6]*/ +10.8170077303506320,
            /*[ 7]*/ -03.0905736372430375,
            /*[ 8]*/ +00.3863217046553797
        );

        //A.ToDebugEnum();
        A.AssertEquals(
            /*[ 0]*/ +01,
            /*[ 1]*/ -06.11292827835928600,
            /*[ 2]*/ +16.52653591409266600,
            /*[ 3]*/ -25.77868616336662600,
            /*[ 4]*/ +25.35080830583167700,
            /*[ 5]*/ -16.08190134129641000,
            /*[ 6]*/ +06.42266125433558000,
            /*[ 7]*/ -01.47559067500681460,
            /*[ 8]*/ +00.14924445948815254
        );

        var filter = new DSP.Filters.ButterworthHighPass(dt, fs, fp, Gp, Gs);

        filter.B.AssertEquals(B);
        filter.A.AssertEquals(A);
    }

    [TestMethod]
    public void TransmissionCoefficient()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ButterworthHighPass(dt, fs, fp, Gp, Gs);

        var transmission__0 = filter.GetTransmissionCoefficient(0, dt);
        var transmission_fp = filter.GetTransmissionCoefficient(fp, dt);
        var transmission_fs = filter.GetTransmissionCoefficient(fs, dt);
        var transmission_fd05 = filter.GetTransmissionCoefficient(fd / 2, dt);

        transmission__0.Abs.AssertLessThan(Gs);
        transmission_fs.Abs.AssertLessOrEqualsThan(Gs);
        transmission_fp.Abs.AssertGreaterOrEqualsThan(Gp, 1.44e-11);
        transmission_fd05.Abs.AssertEquals(1, 1e-15);
    }

    [TestMethod]
    public void SignalProcessing()
    {
        const double fd = 10;                 // Гц // Частота дискретизации
        const double dt = 1 / fd;               // 2с // Период дискретизации

        const double fs = 2 / Consts.pi2;   // Гц // Граничная частота полосы пропускания
        const double fp = 4 / Consts.pi2;   // Гц // Граничная частота полосы запирания

        const double Rp = 1;                    // Неравномерность в полосе пропускания (дБ)
        const double Rs = 40;                   // Неравномерность в полосе пропускания (дБ)

        var Gp = (-Rp).From_dB();
        var Gs = (-Rs).From_dB();

        var filter = new DSP.Filters.ButterworthHighPass(dt, fs, fp, Gp, Gs);

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
        k_s99_db.AssertLessThan(-Rs);               //  \ |    |
        k_s_db.AssertLessOrEqualsThan(-Rs);         //   \+    | fs
                                                    //    \    |
        k_sp_db.AssertGreaterThan(-Rs);             //     \   |
        k_sp_db.AssertLessThan(-Rp);                //      \  |
        k_ps_db.AssertGreaterThan(-Rs);             //       \ |
        k_ps_db.AssertLessThan(-Rp);                //        \+ fp
                                                    //         \0дБ
        k_p_db.AssertEquals(-Rp, 3.6e-3);           //          |
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
        h_sp.AssertGreaterThan(-Rs);                //     \   |
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