﻿using System.Collections;
using System.Diagnostics;

using MathCore.DSP.Filters;

// ReSharper disable RedundantArgumentDefaultValue
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters.ComplexTests;

[TestClassHandler("FailResultHandler")]
public class LowPass_fp500_fs1500_df5000 : ComplexTest
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

    private static class Information
    {
        /// <summary>Частота дискретизации, Гц</summary>
        public const double fd = 5000;

        /// <summary>Период дискретизации</summary>
        public const double dt = 1 / fd;

        /// <summary>Граничная частота полосы пропускания</summary>
        public const double fp = 500;

        /// <summary>Граничная частота полосы запирания</summary>
        public const double fs = 1500;

        /// <summary>Неравномерность в полосе пропускания (дБ)</summary>
        public const double Rp = 1;

        /// <summary>Неравномерность в полосе пропускания (дБ)</summary>
        public const double Rs = 30;

        public static AnalogBasedFilter.Specification GetSpecification() => new(dt, fp, fs, (-Rp).From_dB(), (-Rs).From_dB());
    }

    [TestMethod]
    public void Butterworth()
    {
        var specification = Information.GetSpecification();
        var filter = new DSP.Filters.ButterworthLowPass(specification);
        CheckFilter(filter, 1000);
    }

    [TestMethod]
    public void ChebyshevI()
    {
        var specification = Information.GetSpecification();
        var filter = new DSP.Filters.ChebyshevLowPass(specification, ChebyshevType.I);
        CheckFilter(filter);
    }

    [TestMethod]
    public void ChebyshevII()
    {
        var specification = Information.GetSpecification();
        var filter = new DSP.Filters.ChebyshevLowPass(specification, ChebyshevType.II);

        CheckChebyshevII(filter);
    }

    [TestMethod]
    public void Elliptic()
    {
        var specification = Information.GetSpecification();
        var filter = new DSP.Filters.EllipticLowPass(specification);

        var H = Enumerable.Range(1, 101).Select(i => filter.fd / 100 * i)
           .Select(f => $"{f,5} {GetTransmigrationCoefficient(filter, f, 1000).In_dB_byPower().Round(2)}")
           .ToArray();
            
        CheckFilter(filter);
    }
}