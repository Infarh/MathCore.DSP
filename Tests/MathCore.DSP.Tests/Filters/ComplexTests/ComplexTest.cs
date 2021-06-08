using System;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Fourier;
using MathCore.DSP.Signals;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters.ComplexTests
{
    public abstract class ComplexTest
    {
        protected static void CheckChebyshevII(ChebyshevLowPass filter, int N = 100)
        {
            Assert.That.Value(filter.FilterType).IsEqual(ChebyshevFilter.ChebyshevType.II);

            var specification = filter.Spec;

            var k0 = GetTransmigrationCoefficient(filter, 0, N).In_dB_byPower();
            Assert.That.Value(k0)
               .GreaterOrEqualsThan(-specification.Rp);

            var kp = GetTransmigrationCoefficient(filter, specification.fp, N).In_dB_byPower();
            Assert.That.Value(kp)
               .LessThan(-specification.Rp);

            var ks = GetTransmigrationCoefficient(filter, specification.fs, N).In_dB_byPower();
            Assert.That.Value(ks)
               .LessOrEqualsThan(-specification.Rs);

            for (var (f0, df) = (specification.fs, specification.fd / 2 / 100); f0 <= specification.fd / 2; f0 += df)
            {
                var k = GetTransmigrationCoefficient(filter, f0, N).In_dB_byPower();
                Assert.That.Value(k)
                   .LessThan(-specification.Rs);
            }
        }

        protected static void CheckFilter(AnalogBasedFilter Filter, int N = 1000)
        {
            var specification = Filter.Spec;

            var fp05 = GetTransmigrationCoefficient(Filter, Filter.fp / 2, 10000).In_dB_byPower();
            var fp95 = GetTransmigrationCoefficient(Filter, Filter.fp * 0.95, 10000).In_dB_byPower();
            var fp = GetTransmigrationCoefficient(Filter, Filter.fp, 10000).In_dB_byPower();
            var fps = GetTransmigrationCoefficient(Filter, (Filter.fp + Filter.fs) / 2, 10000).In_dB_byPower();
            var fs = GetTransmigrationCoefficient(Filter, Filter.fs, 10000).In_dB_byPower();
            var fsd = GetTransmigrationCoefficient(Filter, (Filter.fs + Filter.fd) / 2, 10000).In_dB_byPower();

            var Hps = DoubleArrayDSPExtensions.GetTransmissionCoefficient(Filter.A, Filter.B, (Filter.fp + Filter.fs) / 2 * Filter.dt).Power.In_dB_byPower();

            Assert.That.Value(fp05).GreaterThan(-specification.Rp, 0.01);
            Assert.That.Value(fp95).GreaterThan(-specification.Rp);
            Assert.That.Value(fp).GreaterThan(-specification.Rp, 1);

            Assert.That.Value(fps).LessThan(-specification.Rp).GreaterThan(-specification.Rs);
            Assert.That.Value(fs).LessThan(-specification.Rs, 1);
            Assert.That.Value(fsd).LessThan(-specification.Rs);

            var H = Enumerable.Range(0, 101).Select(i => i * Filter.fd / 2 / 100)
               .Select(f => ($"{f,5}", GetTransmigrationCoefficient(Filter, f, 1000).In_dB_byPower()))
               .ToArray();

            CheckBandPass(Filter, specification, N);
            CheckBandStop(Filter, specification, N);
        }

        protected static double GetTransmigrationCoefficient(AnalogBasedFilter Filter, double f0, int N)
        {
            var x = f0 is 0d
                ? MathEnumerableSignal.Const(1, Filter.dt, N)
                : MathEnumerableSignal.Sin(1, f0, 0, Filter.dt, N);

            var y = Filter.ProcessIndividual(x);

            var yy = y.ToArray();

            return y.Power / x.Power;
        }

        /// <summary>Проверка фильтра в полосе пропускания</summary>
        /// <param name="Filter">Проверяемый фильтр</param>
        /// <param name="Specification"></param>
        /// <param name="N"></param>
        private static void CheckBandPass(AnalogBasedFilter Filter, AnalogBasedFilter.Specification Specification, int N)
        {
            var k0 = GetTransmigrationCoefficient(Filter, 0, N).In_dB_byPower();
            Assert.That.Value(k0).GreaterThan(-Specification.Rp);

            var df = Specification.fp / 10;
            for (var f0 = df; Specification.fp - f0 > df; f0 += df)
            {
                var k = GetTransmigrationCoefficient(Filter, f0, N).In_dB_byPower();
                Assert.That.Value(k).GreaterThan(-Specification.Rp, 0.5, $"f0:{f0} fp:{Specification.fp}");
            }

            var kp = GetTransmigrationCoefficient(Filter, Specification.fp, N).In_dB_byPower();

            Assert.That.Value(kp).IsEqual(-Specification.Rp, 0.1);

            var kps_01 = GetTransmigrationCoefficient(Filter, Specification.fp + (Specification.fs - Specification.fp) / 10, N).In_dB_byPower();
            Assert.That.Value(kps_01).LessThan(-Specification.Rp).GreaterThan(-Specification.Rs);
        }

        private static void CheckBandStop(AnalogBasedFilter Filter, AnalogBasedFilter.Specification Specification, int N)
        {
            var kps_05 = GetTransmigrationCoefficient(Filter, Specification.fp + (Specification.fs - Specification.fp) / 2, N).In_dB_byPower();
            Assert.That.Value(kps_05).LessThan(-Specification.Rp).GreaterThan(-Specification.Rs);

            var ks = GetTransmigrationCoefficient(Filter, Specification.fs, N).In_dB_byPower();
            Assert.That.Value(ks).LessOrEqualsThan(-Specification.Rs);

            var df = Specification.fp / 10;
            for (var f0 = Specification.fs; f0 < Specification.fd / 2; f0 += df)
            {
                var k = GetTransmigrationCoefficient(Filter, f0, N).In_dB_byPower();
                Assert.That.Value(k).LessThan(-Specification.Rs);
            }
        }
    }
}