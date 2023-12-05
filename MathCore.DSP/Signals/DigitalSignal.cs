using System.Collections;

using MathCore.DSP.Fourier;
using MathCore.DSP.Signals.Operations;
// ReSharper disable UnusedMember.Global
// ReSharper disable VirtualMemberNeverOverridden.Global

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Signals;

/// <summary>Цифровой сигнал</summary>
[System.Diagnostics.CodeAnalysis.SuppressMessage("Стиль", "IDE1006:Стили именования")]
public abstract class DigitalSignal : IEnumerable<double>, ICollection
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
    public double t0 { get => _t0; init => _t0 = value; }

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

    /// <summary>Мощность - средняя сумма квадратов амплитуд отсчётов</summary>
    public virtual double Power => this.Average(s => s * s);

    /// <summary>Среднее значение отсчётов сигнала - постоянная составляющая</summary>
    public virtual double Average => GetSamples().Average();

    /// <summary>Дисперсия значений сигнала - средняя сумма квадратов отсчётов без квадрата среднего значения</summary>
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
            var (t, dt) = (_t0, _dt);
            foreach (var sample in GetSamples())
                yield return new SignalSample(t += dt, sample);
        }
    }

    /// <summary>Инициализация нового цифрового сигнала</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="t0">Смещение сигнала во времени</param>
    protected DigitalSignal(double dt, double t0 = 0)
    {
        _dt = dt > 0 ? dt : throw new ArgumentOutOfRangeException(nameof(dt), "Период дискретизации должен быть больше 0");
        _t0 = t0;
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
    public virtual EnumerableSignal GetIntegral(double s0 = 0) => new(dt, GetIntegralSamples(s0));

    /// <summary>Получить отсчёты сигнала в виде массива</summary>
    /// <returns>Массив отсчётов сигнала</returns>
    public virtual double[] GetSamples() => this.ToArray();

    public virtual DigitalSpectrum GetSpectrum()
    {
        var spectrum_samples = this.ToArray().FastFourierTransform();
        return new(1 / _dt, spectrum_samples, t0);
    }

    /// <summary>Копирование отсчётов сигнала в массив</summary>
    /// <param name="Destination">Массив места назначения</param>
    /// <param name="Index">Начальный индекс в массиве места назначения</param>
    /// <param name="Length">Длина копируемого участка</param>
    public virtual void CopyTo(double[] Destination, int Index, int Length)
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

    public void CopyTo(double[] Destination, int Index) => CopyTo(Destination, Index, SamplesCount);

    /// <inheritdoc />
    public abstract IEnumerator<double> GetEnumerator();

    /// <inheritdoc />
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    #region Операторы

    public static SumOfSignalsResultSignal operator +(DigitalSignal s1, DigitalSignal s2) => new(s1, s2);
    public static SubstractionOfSignalsResultSignal operator -(DigitalSignal s1, DigitalSignal s2) => new(s1, s2);
    public static MultiplyOfSignalsResultSignal operator *(DigitalSignal s1, DigitalSignal s2) => new(s1, s2);
    public static DivisionOfSignalsResultSignal operator /(DigitalSignal s1, DigitalSignal s2) => new(s1, s2);

    public static SumOfSignalWithScalarResultSignal operator +(DigitalSignal s, double x) => new(s, x);
    public static SumOfSignalWithScalarResultSignal operator +(double x, DigitalSignal s) => new(s, x);
    public static SubstractionOfSignalWithScalarResultSignal operator -(DigitalSignal s, double x) => new(s, x);
    public static SubstractionOfScalarWithSignalResultSignal operator -(double x, DigitalSignal s) => new(s, x);
    public static MultiplyOfSignalWithScalarResultSignal operator *(DigitalSignal s, double x) => new(s, x);
    public static MultiplyOfSignalWithScalarResultSignal operator *(double x, DigitalSignal s) => new(s, x);
    public static DivisionOfSignalWithScalarResultSignal operator /(DigitalSignal s, double x) => new(s, x);
    public static DivisionOfScalarWithSignalResultSignal operator /(double x, DigitalSignal s) => new(s, x);

    public static explicit operator double[](DigitalSignal signal)
    {
        var count = signal.SamplesCount;
        var samples = new double[count];
        signal.CopyTo(samples, 0, count);
        return samples;
    }

    #endregion

    //public virtual Complex[] GetSpectrumSamples(double dt) => FT.FourierTransform()

    #region ICollection

    int ICollection.Count => SamplesCount;

    bool ICollection.IsSynchronized => false;

    object? ICollection.SyncRoot => null;

    void ICollection.CopyTo(Array array, int index)
    {
        if (array is double[] double_array)
        {
            CopyTo(double_array, index);
            return;
        }

        var destination_length = array.Length;
        var i = 0;
        var count = SamplesCount;
        foreach (var sample in this)
        {
            array.SetValue(sample, i);
            i++;
            count--;
            if (i >= destination_length || count == 0) break;
        }
    }

    #endregion
}