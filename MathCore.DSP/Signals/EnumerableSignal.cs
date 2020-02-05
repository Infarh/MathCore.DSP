using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;

namespace MathCore.DSP.Signals
{
    /// <summary>Цифровой сигнал на основе последовательности (потенциально бесконечной) отсчётов</summary>
    public class EnumerableSignal : DigitalSignal
    {
        /// <summary>Перечисление отсчётов сигнала</summary>
        [NotNull] private readonly IEnumerable<double> _Samples;

        /// <summary>Количество отсчётов</summary>
        public override int SamplesCount => _Samples.Count();

        public override double this[int n]
        {
            get => _Samples switch
            {
                double[] array => array[n],
                IList<double> list => list[n],
                _ => _Samples.Skip(n).First()
            };
            set
            {
                switch (_Samples)
                {
                    case double[] array: array[n] = value; break;
                    case IList<double> list: list[n] = value; break;
                    default: throw new NotSupportedException();
                }
            }
        }

        /// <summary>Инициализация нового цифрового сигнал на основе перечисления отсчётов</summary>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="Samples">Перечисление отсчётов сигнала</param>
        public EnumerableSignal(double dt, [NotNull] IEnumerable<double> Samples) : base(dt) => _Samples = Samples ?? throw new ArgumentNullException(nameof(Samples));

        private static IEnumerable<double> GetIntegralSamples(IEnumerable<double> samples, double dt, double s0)
        {
            var dt05 = dt / 2;
            var s = s0;
            yield return s0;
            var last = double.NaN;
            foreach (var sample in samples)
                if (double.IsNaN(last))
                    last = sample;
                else
                    yield return s += (last + (last = sample)) / dt05;
        }

        /// <summary>Вычисление интеграла</summary>
        /// <param name="s0">Константа интегрирования</param>
        /// <returns>Цифровой сигнал, как результат интегрирования</returns>
        public override EnumerableSignal GetIntegral(double s0 = 0) => new EnumerableSignal(_dt, GetIntegralSamples(_Samples, _dt, s0));


        /// <summary>Преобразование в цифровой сигнал на основе массива отсчётов</summary>
        /// <returns>Сигнал на основе массива отсчётов</returns>
        [NotNull]
        public SamplesDigitalSignal ToSamplesSignal() => new SamplesDigitalSignal(_dt, _Samples);

        public override IEnumerator<double> GetEnumerator() => _Samples.GetEnumerator();
    }
}