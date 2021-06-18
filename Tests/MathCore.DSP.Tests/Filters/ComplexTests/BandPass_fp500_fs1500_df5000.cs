using System;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

using Microsoft.VisualStudio.TestTools.UnitTesting;
// ReSharper disable UnusedVariable
// ReSharper disable InconsistentNaming
// ReSharper disable MemberCanBePrivate.Local

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

            var spec = filter.Spec;

            //var hh = Enumerable.Range(1, 151).Select(i => i * 1500d / 150)
            //   .Select(f => (f, K: GetTransmigrationCoefficient(filter, f, 1000)))
            //   .Select(v => $"{v.f,5} {v.K.In_dB_byPower().Round(2)}")
            //   .ToArray();

            CheckTransmissionGreaterThan(filter, -spec.Rp, filter.F, 0.1);

            var h____0 = GetTransmigrationCoefficient(filter, 0, 1000).In_dB_byPower();
            var h___01 = GetTransmigrationCoefficient(filter, 0.1, 10000).In_dB_byPower();
            var h____1 = GetTransmigrationCoefficient(filter, 1, 1000).In_dB_byPower();
            var h_down = GetTransmigrationCoefficient(filter, 10, 1000).In_dB_byPower();
            var h__low = GetTransmigrationCoefficient(filter, 100, 1000).In_dB_byPower();
        }

        public void Elliptic()
        {
        
        }
    }
}