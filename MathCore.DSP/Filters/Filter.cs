using System;
using System.Collections.Generic;
using System.Linq;
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
        public virtual IEnumerable<double> Process(IEnumerable<double> Samples) => Samples.Select(Process);

        /// <summary>Обработать цифровой сигнал</summary>
        /// <param name="Signal">Обрабатываемый сигнал</param>
        /// <returns>Обработанный цифровой сигнал</returns>
        public DigitalSignal Process(DigitalSignal Signal) => new SamplesDigitalSignal((Signal ?? throw new ArgumentNullException(nameof(Signal))).dt, Process((IEnumerable<double>)Signal));

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
        public static DigitalSignal operator *(Filter filter, DigitalSignal signal)
        {
            if (filter is null) throw new ArgumentNullException(nameof(filter));
            return filter.Process(signal ?? throw new ArgumentNullException(nameof(signal)));
        }

        /// <summary>Оператор структурного умножения фильтра и цифрового сигнала</summary>
        /// <param name="filter">Фильтр</param>
        /// <param name="signal">Фильтруемый сигнал</param>
        /// <returns>Сигнал на выходе фильтра</returns>
        public static DigitalSignal operator *(DigitalSignal signal, Filter filter) => filter * signal;

        /// <summary>Оператор структурного параллельного сложения двух фильтров</summary>
        /// <param name="filter1">Первый структурно-объединяемый фильтр</param>
        /// <param name="filter2">Второй структурно-объединяемый фильтр</param>
        /// <returns>Фильтр, представляющий собой параллельное включение двух исходных фильтров</returns>
        public static ParallelFilter operator +(Filter filter1, Filter filter2) => new(filter1, filter2);

        /// <summary>Оператор структурного последовательного сложения двух фильтров</summary>
        /// <param name="filter1">Первый структурно-объединяемый фильтр</param>
        /// <param name="filter2">Второй структурно-объединяемый фильтр</param>
        /// <returns>Фильтр, представляющий собой последовательное включение двух исходных фильтров</returns>
        public static SerialFilter operator *(Filter filter1, Filter filter2) => new(filter1, filter2);
    }
}