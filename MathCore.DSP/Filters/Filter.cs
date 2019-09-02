using System;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;
using MathCore.DSP.Signals;

namespace MathCore.DSP.Filters
{
    /// <summary>Фильтр</summary>
    public abstract class Filter
    {
        /// <summary>Обработать очередной отсчёт сигнала</summary>
        /// <param name="Sample">Обрабатываемый отсчёт сигнала</param>
        /// <returns>Значение на выходе фильтра</returns>
        public abstract double Process(double Sample);

        /// <summary>Обработать последовательность отсчётов сигнала</summary>
        /// <param name="Samples">Перечисление отсчётов сигнала, подлежащая фильтрации</param>
        /// <returns>Последовательность обработанных значений</returns>
        [NotNull] public virtual IEnumerable<double> Process([NotNull] IEnumerable<double> Samples) => Samples.Select(Process);

        /// <summary>Обработать цифровой сигнал</summary>
        /// <param name="Signal">Обрабатываемый сигнал</param>
        /// <returns>Обработанный цифровой сигнал</returns>
        [NotNull] public DigitalSignal Process([NotNull] DigitalSignal Signal) => new SamplesDigitalSignal((Signal ?? throw new ArgumentNullException(nameof(Signal))).dt, Process((IEnumerable<double>)Signal));

        /// <summary>Сбросить состояние фильтра</summary>
        public abstract void Reset();

        /// <summary>Получить комплексный коэффициент передачи</summary>
        /// <param name="f">Частота расчёта комплексного коэффициента передачи</param>
        /// <returns>Значение комплексного коэффициента передачи фильтра на заданной частоте</returns>
        public abstract Complex GetTransmissionCoefficient(double f);

        /// <summary>Оператор структурного умножения фильтра и цифрового сигнала</summary>
        /// <param name="filter">Фильтр</param>
        /// <param name="signal">Фильтруемый сигнал</param>
        /// <returns>Сигнал на выходе фильтра</returns>
        [NotNull]
        public static DigitalSignal operator *([NotNull] Filter filter, [NotNull] DigitalSignal signal)
        {
            if (filter is null) throw new ArgumentNullException(nameof(filter));
            if (signal is null) throw new ArgumentNullException(nameof(signal));
            return filter.Process(signal);
        }

        /// <summary>Оператор структурного умножения фильтра и цифрового сигнала</summary>
        /// <param name="filter">Фильтр</param>
        /// <param name="signal">Фильтруемый сигнал</param>
        /// <returns>Сигнал на выходе фильтра</returns>
        [NotNull] public static DigitalSignal operator *([NotNull] DigitalSignal signal, [NotNull] Filter filter) => filter * signal;

        /// <summary>Оператор структурного параллельного сложения двух фильтров</summary>
        /// <param name="filter1">Первый структурно-объединяемый фильтр</param>
        /// <param name="filter2">Второй структурно-объединяемый фильтр</param>
        /// <returns>Фильтр, представляющий собой параллельное включение двух исходных фильтров</returns>
        [NotNull] public static ParallelFilter operator +([NotNull] Filter filter1, [NotNull] Filter filter2) => new ParallelFilter(filter1, filter2);

        /// <summary>Оператор структурного последовательного сложения двух фильтров</summary>
        /// <param name="filter1">Первый структурно-объединяемый фильтр</param>
        /// <param name="filter2">Второй структурно-объединяемый фильтр</param>
        /// <returns>Фильтр, представляющий собой последовательное включение двух исходных фильтров</returns>
        [NotNull] public static SerialFilter operator *([NotNull] Filter filter1, [NotNull] Filter filter2) => new SerialFilter(filter1, filter2);
    }
}