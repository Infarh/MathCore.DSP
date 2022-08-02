using System;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;

// ReSharper disable RedundantArgumentDefaultValue

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters.ComplexTests;

[TestClass]
public class NormedLowPass : ComplexTest
{
    private static class Information
    {
        /// <summary>Частота дискретизации, Гц</summary>
        public const double fd = 5;

        /// <summary>Период дискретизации</summary>
        public const double dt = 1 / fd;

        /// <summary>Граничная частота полосы пропускания</summary>
        public const double fp = 1;

        /// <summary>Граничная частота полосы запирания</summary>
        public const double fs = 2;

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
        CheckFilter(filter);
    }

    [TestMethod]
    public void ChebyshevI()
    {
        var specification = Information.GetSpecification();
        var filter = new ChebyshevLowPass(specification, ChebyshevType.I);
        CheckFilter(filter);
    }

    [TestMethod]
    public void ChebyshevII()
    {
        var specification = Information.GetSpecification();
        var filter = new ChebyshevLowPass(specification, ChebyshevType.II);

        CheckChebyshevII(filter);
    }

    [TestMethod]
    public void Elliptic()
    {
        var specification = Information.GetSpecification();
        var filter = new DSP.Filters.EllipticLowPass(specification);
        CheckFilter(filter);
    }
}