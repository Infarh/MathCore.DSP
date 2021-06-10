using System;
using System.Linq;

using MathCore.DSP.Filters;
using MathCore.DSP.Signals;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters.ComplexTests
{
    /// <summary>Базовый интеграционный тест цифрового фильтра</summary>
    public abstract class ComplexTest
    {
        /// <summary>Тест фильтра Чебышева II-рода</summary>
        /// <param name="filter">Тестируемый фильтр чебышева</param>
        /// <param name="N">Число точек сигнала</param>
        protected static void CheckChebyshevII(ChebyshevLowPass filter, int N = 100)
        {
            // Фильтр должен быть фильтром Чебышева второго рода
            Assert.That.Value(filter.FilterType).IsEqual(ChebyshevFilter.ChebyshevType.II);

            var specification = filter.Spec;

            // Коэффициент передачи на нулевой частоте должен быть не больше коэффициента подавления в полосе пропускания
            var k0 = GetTransmigrationCoefficient(filter, 0, N).In_dB_byPower();
            Assert.That.Value(k0).GreaterOrEqualsThan(-specification.Rp);

            // Коэффициент передачи на граничной частоте полосы пропускания должно быть меньше коэффициента передачи в полосе пропускания
            var kp = GetTransmigrationCoefficient(filter, specification.fp, N).In_dB_byPower();
            Assert.That.Value(kp).LessThan(-specification.Rp);

            // Коэффициент передачи на граничной частоте полосы заграждения должен быть меньше коэффициента пропускания полосы подавления
            var ks = GetTransmigrationCoefficient(filter, specification.fs, N).In_dB_byPower();
            Assert.That.Value(ks).LessOrEqualsThan(-specification.Rs);

            // Коэффициенты передачи в полосе заграждения должен быть не больше чем коэффициент передачи полосы заграждения
            for (var (f0, df) = (specification.fs, specification.fd / 2 / 100); f0 <= specification.fd / 2; f0 += df)
            {
                var k = GetTransmigrationCoefficient(filter, f0, N).In_dB_byPower();
                Assert.That.Value(k).LessThan(-specification.Rs);
            }
        }


        private static void CheckTransmissionGreaterThan(AnalogBasedFilter Filter, double Gdb, double f1, double f2, double Accuracy = 0, int M = 100, int N = 10000)
        {
            var df = (f2 - f1) / (M + 1);
            for (var i = 0; i < M; i++)
            {
                var f = f1 + i * df;
                var Hdb = GetTransmigrationCoefficient(Filter, f, N).In_dB_byPower();
                if (Accuracy is 0d)
                    Assert.That.Value(Hdb).GreaterThan(Gdb, $"Коэффициент передачи на частоте {f} должен быть меньше {Gdb}, а был {Hdb}, при точности {Accuracy}");
                else
                    Assert.That.Value(Hdb).GreaterThan(Gdb, Accuracy, $"Коэффициент передачи на частоте {f} должен быть меньше {Gdb}, а был {Hdb}");
            }
        }

        private static void CheckTransmissionLessThen(AnalogBasedFilter Filter, double Gdb, double f1, double f2, double Accuracy = 0, int M = 100, int N = 10000)
        {
            var df = (f2 - f1) / (M + 1);
            for (var i = 0; i < M; i++)
            {
                var f = f1 + i * df;
                var Hdb = GetTransmigrationCoefficient(Filter, f, N).In_dB_byPower();
                if (Accuracy is 0d)
                    Assert.That.Value(Hdb).LessOrEqualsThan(Gdb, $"Коэффициент передачи на частоте {f} должен быть меньше {Gdb}, а был {Hdb}, при точности {Accuracy}");
                else
                    Assert.That.Value(Hdb).LessOrEqualsThan(Gdb, Accuracy, $"Коэффициент передачи на частоте {f} должен быть меньше {Gdb}, а был {Hdb}");
            }
        }

        protected static void CheckFilter(AnalogBasedFilter Filter, int N = 1000)
        {
            var specification = Filter.Spec;

            var fp05 = GetTransmigrationCoefficient(Filter, specification.fp / 2, 10000).In_dB_byPower();
            var fp95 = GetTransmigrationCoefficient(Filter, specification.fp * 0.95, 10000).In_dB_byPower();
            var fp = GetTransmigrationCoefficient(Filter, specification.fp, 10000).In_dB_byPower();
            var fps = GetTransmigrationCoefficient(Filter, (specification.fp + specification.fs) / 2, 10000).In_dB_byPower();
            var fs = GetTransmigrationCoefficient(Filter, specification.fs, 10000).In_dB_byPower();
            var fsd = GetTransmigrationCoefficient(Filter, (specification.fs + specification.fd / 2) / 2, 10000).In_dB_byPower();

            CheckTransmissionGreaterThan(Filter, -specification.Rp, 0, specification.fp, 0.1);
            CheckTransmissionLessThen(Filter, -specification.Rp, specification.fp, specification.fs);
            CheckTransmissionLessThen(Filter, -specification.Rp, specification.fs, specification.fd / 2);

            //var Hps = DoubleArrayDSPExtensions.GetTransmissionCoefficient(Filter.A, Filter.B, (Filter.fp + Filter.fs) / 2 * Filter.dt).Power.In_dB_byPower();

            Assert.That.Value(fp05).GreaterThan(-specification.Rp);       // В середине полосы пропускания коэффициент передачи не должен быть ниже Rp
            Assert.That.Value(fp95).GreaterThan(-specification.Rp);       // В конце полосы пропускания коэффициент передачи не должен быть ниже Rp
            Assert.That.Value(fp).IsEqual(-specification.Rp, 4.8e-3);     // На границе полосы пропускания

            Assert.That.Value(fps).LessThan(-specification.Rp);           // Между граничными частотам полосы пропускания и полосы подавления коэффициент передачи фильтра должен быть ниже Rp
            Assert.That.Value(fs).LessThan(-specification.Rs, 1);         // На граничной частоте полосы подавления коэффициент передачи должен быть ниже Rp
            Assert.That.Value(fsd).LessThan(-specification.Rs);           // В полосе подавления (между граничной частотой полосы подавления и половиной частоты дискретизации) коэффициент передачи должен быть ниже Rp

            CheckBandPass(Filter, specification, N);
            CheckBandStop(Filter, specification, N);
        }

        protected static double GetTransmigrationCoefficient(AnalogBasedFilter Filter, double f0, int N)
        {
            var x = f0 is 0d
                ? MathEnumerableSignal.Const(1, Filter.dt, N)
                : MathEnumerableSignal.Sin(1, f0, 0, Filter.dt, N);

            return Filter.ProcessIndividual(x).Power / x.Power;
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
            Assert.That.Value(kps_05).LessThan(-Specification.Rp);

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