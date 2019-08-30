using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;
using MathCore.DSP.Signals;

namespace MathCore.DSP.Filters
{
    /// <summary>Цифровой фильтр</summary>
    public abstract class DigitalFilter : Filter
    {
        /// <summary>Преобразование частоты цифрового фильтра в частоту аналогового прототипа</summary>
        /// <param name="DigitalFrequency">Значение на оси частот цифрового фильтра</param>
        /// <param name="dt">Период дискретизации</param>
        /// <returns>Значение на оси частот аналогового прототипа</returns>
        public static double ToAnalogFrequency(double DigitalFrequency, double dt) => Math.Tan(Math.PI * DigitalFrequency * dt) / (Math.PI * dt);

        public static Complex ToZ(Complex p, double dt) => (2 / dt + p) / (2 / dt - p);

        public static double GetNomalizeCoefficient([NotNull] IEnumerable<Complex> poles, double dt)
        {
            if (poles is null) throw new ArgumentNullException(nameof(poles));

            var result = new Complex(1);
            var k = 2 / dt;
            result = poles.Aggregate(result, (current, p0) => current * (k - p0));
            if (Math.Abs(result.Im / result.Re) > 1e-15) throw new InvalidOperationException("Комплексный результат");
            return 1 / result.Re;
        }

        /// <summary>Вектор состояния</summary>
        protected readonly double[] State;

        /// <summary>Порядок фильтра</summary>
        public virtual int Order => State.Length;

        /// <summary>Инициализация нового цифрового фильтра</summary>
        /// <param name="Order">Порядок фильтра</param>
        protected DigitalFilter(int Order) => State = new double[Order];

        /// <summary>Обработать очередной отсчёт цифрового сигнала</summary>
        /// <param name="sample">Обрабатываемый отсчёт цифрового сигнала</param>
        /// <param name="state">Вектор состояния фильтра</param>
        /// <returns>Значение сигнала на выходе фильтра после обработки отсчёта</returns>
        public abstract double Process(double sample, [NotNull] double[] state);

        /// <summary>Обработать отсчёт цифрового сигнала</summary>
        /// <param name="sample">Обрабатываемый отсчёт цифрового сигнала</param>
        /// <returns>Значение сигнала на выходе фильтра после обработки отсчёта</returns>
        public override double Process(double sample) => Process(sample, State);

        [NotNull]
        public DigitalSignal Process([NotNull] DigitalSignal Signal, [NotNull] double[] state)
        {
            if (Signal is null) throw new ArgumentNullException(nameof(Signal));
            if (state is null) throw new ArgumentNullException(nameof(state));
            if (state.Length != Order) throw new InvalidOperationException($"Длина вектора состояний {state.Length} не равна порядку фильтра {Order}");

            return new SamplesDigitalSignal(Signal.dt, Signal.Select(s => Process(s, state)));
        }

        [NotNull] public DigitalSignal ProcessIndividual([NotNull] DigitalSignal Signal) => Process(Signal, new double[Order]);

        /// <summary>Сбросить состояние фильтра</summary>
        public override void Reset() => Array.Clear(State, 0, State.Length);

        public Complex GetTransmissionCoefficient(double f, double dt) => GetTransmissionCoefficient(f * dt);
    }
}
