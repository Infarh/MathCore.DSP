using System;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters.ComplexTests
{
    [TestClass]
    public class BandPass_fp500_fs1500_df5000 : ComplexTest
    {
        private static class Information
        {
            /// <summary>Частота дискретизации, Гц</summary>
            public const double fd = 5000;

            /// <summary>Период дискретизации</summary>
            public const double dt = 1 / fd;

            /// <summary>Граничная частота полосы пропускания</summary>
            public const double fp = 10;

            /// <summary>Граничная частота полосы запирания</summary>
            public const double fs = 40;

            /// <summary>Неравномерность в полосе пропускания (дБ)</summary>
            public const double Rp = 1;

            /// <summary>Неравномерность в полосе пропускания (дБ)</summary>
            public const double Rs = 30;

            /// <summary>Нижняя граница полосы пропускания</summary>
            public const double fmin = 100;

            /// <summary>Верхняя граница полосы пропускания</summary>
            public const double fmax = 1000;

            public static AnalogBasedFilter.Specification GetSpecification() => new(dt, fp, fs, (-Rp).From_dB(), (-Rs).From_dB());
        }

        [TestMethod]
        public void Butterworth()
        {
            var filter = new ButterworthBandPass(
                Information.dt, 
                Information.fp, Information.fs, 
                Information.fmin, Information.fmax, 
                (-Information.Rp).From_dB(), (-Information.Rs).From_dB());

        }
    }
}