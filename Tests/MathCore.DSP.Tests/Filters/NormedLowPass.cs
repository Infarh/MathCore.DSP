using System;
using System.Collections.Generic;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;
using MathCore.DSP.Signals.Implementations.Enumerables.Base;
using MathCore.Vectors;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using ButterworthLowPass = MathCore.DSP.Filters.ButterworthLowPass;
// ReSharper disable RedundantArgumentDefaultValue

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Tests.Filters
{
    [TestClass]
    public class NormedLowPass
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
            var filter = new DSP.Filters.ButterworthLowPass(Information.GetSpecification());
            CheckFilter(filter);
        }

        [TestMethod]
        public void ChebyshevI()
        {
            var filter = new ChebyshevLowPass(Information.GetSpecification(), ChebyshevFilter.ChebyshevType.I);
            CheckFilter(filter);
        }

        [TestMethod]
        public void ChebyshevII()
        {
            var specification = Information.GetSpecification();
            var filter = new ChebyshevLowPass(specification, ChebyshevFilter.ChebyshevType.II);

            var k = new List<Vector2D>();
            for (var (f0, df) = (0d, Information.fd / 2 / 100); f0 <= Information.fd / 2; f0 += df)
                k.Add(new(f0, GetTransmitionCoefficient(filter, f0).In_dB_byPower()));

            //GetTransmitionCoefficient(filter, specification.fs).In_dB_byPower()
        }

        [TestMethod]
        public void Elliptic()
        {
            var filter = new DSP.Filters.EllipticLowPass(Information.GetSpecification());
            CheckFilter(filter);
        }

        private void CheckFilter(AnalogBasedFilter Filter)
        {
            //var k = new List<double>();
            //for (var (f0, df) = (0d, Information.fd / 2 / 100); f0 <= Information.fd / 2; f0 += df)
            //    k.Add(GetTransmitionCoefficient(Filter, f0));

            CheckBandPass(Filter);
            CheckBandStop(Filter);
        }

        private static double GetTransmitionCoefficient(DigitalFilter Filter, double f0, int N = 100)
        {
            var x = f0 is 0d
                ? MathEnumerableSignal.Const(1, Information.dt, N)
                : MathEnumerableSignal.Sin(1, f0, 0, Information.dt, N);
            var y = Filter.ProcessIndividual(x);
            return y.Power / x.Power;
        }

        private static void CheckBandPass(DigitalFilter Filter)
        {
            const double df = Information.fp / 10;
            var k0 = GetTransmitionCoefficient(Filter, 0).In_dB_byPower();
            Assert.That.Value(k0)
               .GreaterThan(-Information.Rp);

            for (var f0 = df; Information.fp - f0 > df; f0 += df)
            {
                var k = GetTransmitionCoefficient(Filter, f0).In_dB_byPower();
                Assert.That.Value(k)
                   .GreaterThan(-3)
                   .GreaterThan(-Information.Rp - 0.1);
            }

            var kp = GetTransmitionCoefficient(Filter, Information.fp).In_dB_byPower();
            Assert.That.Value(kp)
               .LessThan(-Information.Rp)
               .GreaterOrEqual(-3);

            var kp1 = GetTransmitionCoefficient(Filter, Information.fp + df).In_dB_byPower();
            Assert.That.Value(kp1)
               .LessThan(-3 + 1);

            var kp2 = GetTransmitionCoefficient(Filter, Information.fp + 2 * df).In_dB_byPower();
            Assert.That.Value(kp2)
               .LessThan(-3);
        }

        private static void CheckBandStop(DigitalFilter Filter)
        {
            const double df = Information.fp / 10;

            var ks = GetTransmitionCoefficient(Filter, Information.fs).In_dB_byPower();
            Assert.That.Value(ks)
               .LessOrEqualsThan(-Information.Rs + 1);

            for (var f0 = Information.fs; f0 < Information.fd / 2; f0 += df)
            {
                var k = GetTransmitionCoefficient(Filter, f0).In_dB_byPower();
                Assert.That.Value(k)
                   .LessThan(-3)
                   .LessThan(-Information.Rs + 1);
            }
        }
    }
}
