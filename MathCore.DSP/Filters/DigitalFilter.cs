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

        /// <summary>Преобразование частоты аналогового  прототипа в частоту цифрового фильтра</summary>
        /// <param name="AnalogFrequency">Значение на оси частот аналогового фильтра</param>
        /// <param name="dt">Период дискретизации</param>
        /// <returns>Значение на оси частот цифрового фильтра</returns>
        public static double ToDigitalFrequency(double AnalogFrequency, double dt) => Math.Atan(Math.PI * AnalogFrequency * dt) / (Math.PI * dt);

        /// <summary>Преобразование полюса из p-плоскости в z-плоскость</summary>
        /// <param name="p">Полюс p-плоскости</param>
        /// <param name="dt">Период дискретизации</param>
        /// <returns>Полюс в z-плоскости</returns>
        public static Complex ToZ(Complex p, double dt)
        {
            var w = 2 / dt;
            return (w + p) / (w - p);
        }

        public static (Complex[] zPoles, Complex[] zZeros) ToZ(Complex[] pPoles, Complex[] pZeros, double dt)
        {
            if (pPoles is null) throw new ArgumentNullException(nameof(pPoles));
            if (pZeros is null) throw new ArgumentNullException(nameof(pZeros));
            if (pZeros.Length > pPoles.Length) throw new ArgumentException("Число нулей не должно превышать числа полюсов", nameof(pZeros));

            var poles_count = pPoles.Length;

            var zZeros = new Complex[poles_count];
            var zPoles = new Complex[poles_count];

            for (var i = 0; i < pPoles.Length; i++) zPoles[i] = ToZ(pPoles[i], dt);
            for (var i = 0; i < pZeros.Length; i++) zZeros[i] = ToZ(pZeros[i], dt);
            for (var i = pZeros.Length; i < zZeros.Length; i++) zZeros[i] = -1;

            return (zPoles, zZeros);
        }

        public static IEnumerable<Complex> ToZ(IEnumerable<Complex> p, double dt) => p.Select(z => ToZ(z, dt));

        public static IEnumerable<Complex> ToZ(IEnumerable<Complex> p, double W0, double dt) => p.Select(z => ToZ(z * W0, dt));

        public static Complex[] ToZArray(IEnumerable<Complex> p, double dt, double W0 = 1) => ToZ(p, W0, dt).ToArray();

        /// <summary>Преобразование полюса из z-плоскости в p-плоскость</summary>
        /// <param name="z">Полюс z-плоскости</param>
        /// <param name="dt">Период дискретизации</param>
        /// <returns>Полюс в p-плоскости</returns>
        public static Complex ToP(Complex z, double dt) => 2 / dt * (z - 1) / (z + 1);

        /// <summary>Расчёт нормирующего множителя (приводящего системную-передаточную функцию к виду с максимумом в 1)</summary>
        /// <param name="poles">Набор полюсов</param>
        /// <param name="dt">Период дискретизации</param>
        /// <returns>Нормирующий множитель</returns>
        public static double GetNormalizeCoefficient([NotNull] IEnumerable<Complex> poles, double dt)
        {
            if (poles is null) throw new ArgumentNullException(nameof(poles));

            var k = 2 / dt;
            var (re, im) = poles.Aggregate(Complex.Real, (current, p0) => current * (k - p0));
            return (im / re).Abs() <= 1e-15
                ? 1 / re
                : throw new InvalidOperationException("Комплексный результат");
        }

        /// <summary>Вектор состояния</summary>
        protected readonly double[] State;

        /// <summary>Порядок фильтра</summary>
        public virtual int Order => State.Length - 1;

        /// <summary>Инициализация нового цифрового фильтра</summary>
        /// <param name="Order">Порядок фильтра</param>
        protected DigitalFilter(int Order) => State = new double[Order];

        /// <summary>Обработать очередной отсчёт цифрового сигнала</summary>
        /// <param name="Sample">Обрабатываемый отсчёт цифрового сигнала</param>
        /// <param name="state">Вектор состояния фильтра</param>
        /// <returns>Значение сигнала на выходе фильтра после обработки отсчёта</returns>
        public abstract double Process(double Sample, [NotNull] double[] state);

        /// <summary>Обработать отсчёт цифрового сигнала</summary>
        /// <param name="Sample">Обрабатываемый отсчёт цифрового сигнала</param>
        /// <returns>Значение сигнала на выходе фильтра после обработки отсчёта</returns>
        public override double Process(double Sample) => Process(Sample, State);

        /// <summary>Обработать цифровой сигнал</summary>
        /// <param name="Signal">Цифровой сигнал</param>
        /// <param name="state">Вектор состояния фильтра</param>
        /// <returns>Обработанный цифровой сигнал</returns>
        [NotNull]
        public DigitalSignal Process([NotNull] DigitalSignal Signal, [NotNull] double[] state)
        {
            if (Signal is null) throw new ArgumentNullException(nameof(Signal));
            if (state is null) throw new ArgumentNullException(nameof(state));
            if (state.Length != Order + 1) throw new InvalidOperationException($"Длина вектора состояний {state.Length} не равна порядку фильтра {Order} + 1");

            return new SamplesDigitalSignal(Signal.dt, Signal.Select(s => Process(s, state)));
        }

        /// <summary>Обработать цифровой сигнал независимо от состояния фильтра (вектор состояния создаётся на каждый вызов этого метода)</summary>
        /// <param name="Signal">Обрабатываемый цифровой сигнал</param>
        /// <returns>Обработанный цифровой сигнал</returns>
        [NotNull] public DigitalSignal ProcessIndividual([NotNull] DigitalSignal Signal) => Process(Signal, new double[Order + 1]);

        /// <summary>Сбросить состояние фильтра</summary>
        public override void Reset() => Array.Clear(State, 0, State.Length);

        /// <summary>Получить коэффициент передачи фильтра на указанной частоте (КЧХ)</summary>
        /// <param name="f">Частота расчёта коэффициента передачи</param>
        /// <param name="dt">Период дискретизации</param>
        /// <returns>Комплексный коэффициент передачи фильтра</returns>
        public Complex GetTransmissionCoefficient(double f, double dt) => GetTransmissionCoefficient(f * dt);
    }
}
