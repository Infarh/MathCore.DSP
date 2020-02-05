using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MathCore.Annotations;
using MathCore.DSP.Signals.Operations;
// ReSharper disable UnusedMember.Global
// ReSharper disable VirtualMemberNeverOverridden.Global

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Signals
{
    /// <summary>Цифровой сигнал</summary>
    [System.Diagnostics.CodeAnalysis.SuppressMessage("Стиль", "IDE1006:Стили именования")]
    public abstract class DigitalSignal : IEnumerable<double>
    {
        /// <summary>Период дискретизации</summary>
        protected readonly double _dt;

        /// <summary>Смещение начала сигнала во времени</summary>
        protected readonly double _t0;

        /// <summary>Период дискретизации</summary>
        public double dt => _dt;

        /// <summary>Частота дискретизации</summary>
        public double fd => 1 / _dt;

        /// <summary>Начальное смещение сигнала во времени</summary>
        public double t0 => _t0;

        /// <summary>Полное время сигнала</summary>
        public virtual double TotalTime => SamplesCount * _dt;

        /// <summary>Количество отсчётов</summary>
        public abstract int SamplesCount { get; }

        /// <summary>Минимальное значение</summary>
        public virtual double Min => this.Min();

        /// <summary>Максимальное значение</summary>
        public virtual double Max => this.Max();

        /// <summary>Амплитуда</summary>
        public virtual double PeakToPeakAmplitude => this.GetMinMax().Length;

        /// <summary>Мощность</summary>
        public virtual double Power => this.Average(s => s * s);

        /// <summary>Среднее значение отсчётов сигнала</summary>
        public virtual double Average => GetSamples().Average();

        /// <summary>Дисперсия значений сигнала</summary>
        public virtual double Variance => GetSamples().Dispersion();

        /// <summary>Отсчёты по индексу</summary>
        /// <param name="n">Индекс отсчёта</param>
        /// <returns>Значение отсчёта</returns>
        public abstract double this[int n] { get; set; }

        /// <summary>Отсчёты сигнала</summary>
        public virtual IEnumerable<SignalSample> Samples
        {
            get
            {
                var t = _t0;
                foreach (var sample in GetSamples())
                    yield return new SignalSample(t += _dt, sample);
            }
        }

        /// <summary>Инициализация нового цифрового сигнала</summary>
        /// <param name="dt">Период дискретизации</param>
        protected DigitalSignal(double dt)
        {
            if (dt <= 0) throw new ArgumentOutOfRangeException(nameof(dt), "Период дискретизации должен быть больше 0");
            _dt = dt;
        }

        /// <summary>Перечисление отсчётов интеграла сигнала</summary>
        /// <param name="s0">Константа интегрирования</param>
        /// <returns>Перечисление отсчётов интеграла сигнала</returns>
        protected virtual IEnumerable<double> GetIntegralSamples(double s0)
        {
            var s = s0;
            yield return s0;
            var dt05 = dt / 2;
            var last = this[0];
            for (int i = 1, count = SamplesCount; i < count; i++)
                yield return s += (last + (last = this[i])) / dt05;
        }

        /// <summary>Вычисление интеграла</summary>
        /// <param name="s0">Константа интегрирования</param>
        /// <returns>Цифровой сигнал, как результат интегрирования</returns>
        [NotNull] public virtual EnumerableSignal GetIntegral(double s0 = 0) => new EnumerableSignal(dt, GetIntegralSamples(s0));

        /// <summary>Получить отсчёты сигнала в виде массива</summary>
        /// <returns>Массив отсчётов сигнала</returns>
        [NotNull] public virtual double[] GetSamples() => this.ToArray();

        /// <summary>Копирование отсчётов сигнала в массив</summary>
        /// <param name="Destination">Массив места назначения</param>
        /// <param name="Index">Начальный индекс в массиве места назначения</param>
        /// <param name="Length">Длина копируемого участка</param>
        public virtual void CopyTo([NotNull] double[] Destination, int Index, int Length)
        {
            var destination_length = Destination.Length;
            if (Index >= destination_length || Length < 1) return;
            var i = Index;
            var count = Length;
            foreach (var sample in this)
            {
                Destination[i] = sample;
                i++;
                count--;
                if (i >= destination_length || count == 0) break;
            }
        }

        /// <inheritdoc />
        public abstract IEnumerator<double> GetEnumerator();

        /// <inheritdoc />
        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

        #region Операторы

        [NotNull] public static SumOfSignalsResultSignal operator +([NotNull] DigitalSignal s1, [NotNull] DigitalSignal s2) => new SumOfSignalsResultSignal(s1, s2);
        [NotNull] public static SubstractionOfSignalsResultSignal operator -([NotNull] DigitalSignal s1, [NotNull] DigitalSignal s2) => new SubstractionOfSignalsResultSignal(s1, s2);
        [NotNull] public static MultiplyOfSignalsResultSignal operator *([NotNull] DigitalSignal s1, [NotNull] DigitalSignal s2) => new MultiplyOfSignalsResultSignal(s1, s2);
        [NotNull] public static DivisionOfSignalsResultSignal operator /([NotNull] DigitalSignal s1, [NotNull] DigitalSignal s2) => new DivisionOfSignalsResultSignal(s1, s2);

        [NotNull] public static SumOfSignalWithScalarResultSignal operator +([NotNull] DigitalSignal s, double x) => new SumOfSignalWithScalarResultSignal(s, x);
        [NotNull] public static SumOfSignalWithScalarResultSignal operator +(double x, [NotNull] DigitalSignal s) => new SumOfSignalWithScalarResultSignal(s, x);
        [NotNull] public static SubstractionOfSignalWithScalarResultSignal operator -([NotNull] DigitalSignal s, double x) => new SubstractionOfSignalWithScalarResultSignal(s, x);
        [NotNull] public static SubstractionOfScalarWithSignalResultSignal operator -(double x, [NotNull] DigitalSignal s) => new SubstractionOfScalarWithSignalResultSignal(s, x);
        [NotNull] public static MultiplyOfSignalWithScalarResultSignal operator *([NotNull] DigitalSignal s, double x) => new MultiplyOfSignalWithScalarResultSignal(s, x);
        [NotNull] public static MultiplyOfSignalWithScalarResultSignal operator *(double x, [NotNull] DigitalSignal s) => new MultiplyOfSignalWithScalarResultSignal(s, x);
        [NotNull] public static DivisionOfSignalWithScalarResultSignal operator /([NotNull] DigitalSignal s, double x) => new DivisionOfSignalWithScalarResultSignal(s, x);
        [NotNull] public static DivisionOfScalarWithSignalResultSignal operator /(double x, [NotNull] DigitalSignal s) => new DivisionOfScalarWithSignalResultSignal(s, x);

        #endregion

        //public virtual Complex[] GetSpectrumSamples(double dt) => FT.FourierTransform()
    }

}