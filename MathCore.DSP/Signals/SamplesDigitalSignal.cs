using MathCore.Annotations;

using Suppress = System.Diagnostics.CodeAnalysis.SuppressMessageAttribute;
// ReSharper disable InconsistentNaming
// ReSharper disable UnusedMember.Global
// ReSharper disable ConvertToAutoPropertyWhenPossible
// ReSharper disable LoopCanBeConvertedToQuery
// ReSharper disable ForCanBeConvertedToForeach

namespace MathCore.DSP.Signals;

/// <summary>Цифровой сигнал на основе массива отсчётов</summary>
[Suppress("Performance", "CA1819:Properties should not return arrays")]
public class SamplesDigitalSignal : DigitalSignal, IEquatable<SamplesDigitalSignal>, IIndexableRef<double>
{
    /// <summary>Массив отсчётов</summary>
    private readonly double[] _Samples;

    /// <inheritdoc />
    /// <summary>Число элементов массива</summary>
    public override int SamplesCount => _Samples.Length;

    /// <inheritdoc />
    /// <summary>Возвращает, или записывает значение элемента массива с указанным индексом</summary>
    public override double this[int n]
    {
        get => _Samples[n];
        set => _Samples[n] = value;
    }

    /// <inheritdoc />
    ref double IIndexableRef<double>.this[int n] => ref _Samples[n];

    /// <inheritdoc />
    public override IEnumerable<SignalSample> Samples
    {
        get
        {
            var (t0, dt, samples, count) = (_t0, _dt, _Samples, _Samples.Length);
            for (var i = 0; i < count; i++)
                yield return new SignalSample(i * dt + t0, samples[i]);
        }
    }

    public SamplesDigitalSignal(double dt, double[] Samples, double t0 = 0) : base(dt, t0) => _Samples = Samples.NotNull();

    public SamplesDigitalSignal(double dt, int SamplesCount, Func<double, double> f, double t0 = 0) : this(dt, f.NotNull().Sampling(0, dt, SamplesCount), t0) { }

    private static double[] GetSamplesArray(IEnumerable<double> Samples) => Samples switch
    {
        double[] array => array,
        { } => Samples.ToArray(),
        _ => throw new ArgumentNullException(nameof(Samples))
    };
    public SamplesDigitalSignal(double dt, IEnumerable<double> Samples, double t0 = 0) : this(dt, GetSamplesArray(Samples.NotNull()), t0) { }

    public SamplesDigitalSignal(double dt, IEnumerable<int> Samples, double t0 = 0) : this(dt, Samples.NotNull().ToArray(v => (double)v), t0) { }

    /// <inheritdoc />
    protected override IEnumerable<double> GetIntegralSamples(double s0)
    {
        var dt05 = dt / 2;
        var s = s0;
        yield return s0;

        var samples = _Samples;
        var last = samples[0];
        for (int i = 1, count = samples.Length; i < count; i++)
            yield return s += (last + (last = samples[i])) / dt05;
    }

    /// <summary>Вычисление интеграла</summary>
    /// <param name="s0">Константа интегрирования</param>
    /// <returns>Цифровой сигнал, как результат интегрирования</returns>
    public SamplesDigitalSignal GetIntegralSampled(double s0 = 0)
    {
        var source = _Samples;
        var source_length = source.Length;
        var samples = new double[source_length];
        var dt05 = _dt / 2;
        var s = s0;
        var last = source[0];
        samples[0] = s0;
        for (var i = 1; i < source_length; i++)
            samples[i] = s += (last + (last = source[i])) / dt05;

        return new SamplesDigitalSignal(_dt, samples);
    }

    /// <inheritdoc />
    public override double[] GetSamples() => _Samples;

    /// <inheritdoc />
    public override void CopyTo(double[] Destination, int Index, int Length) => Array.Copy(_Samples, 0, Destination, Index, Length);

    #region Override Object

    /// <inheritdoc />
    public override string ToString() => $"signal dt:{_dt}; count:{_Samples.Length}; power:{Power.RoundAdaptive(2)}";

    /// <inheritdoc />
    public bool Equals(SamplesDigitalSignal? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;

        var samples = _Samples;
        var length = samples.Length;
        if (!_dt.Equals(other._dt) || length != other._Samples.Length) return false;
        if (ReferenceEquals(samples, other._Samples)) return true;

        for (var i = 0; i < length; i++)
            if (samples[i] != other._Samples[i])
                return false;

        return true;
    }

    /// <inheritdoc />
    public override bool Equals(object? obj) =>
        obj != null
        && (ReferenceEquals(this, obj)
            || obj.GetType() == GetType() && Equals((SamplesDigitalSignal)obj));

    /// <inheritdoc />
    public override int GetHashCode()
    {
        var hash = _dt.GetHashCode() * 397;
        hash = unchecked((hash * 397) ^ _t0.GetHashCode());
        var samples = _Samples;
        var length = samples.Length;
        for (var i = 0; i < length; i++)
            hash = unchecked((hash * 397) ^ samples[i].GetHashCode());
        return hash;
    }

    #endregion

    /// <inheritdoc />
    public override IEnumerator<double> GetEnumerator() => ((IEnumerable<double>)_Samples).GetEnumerator();

    /// <summary>Массив неявного преобразования типа <see cref="SamplesDigitalSignal"/> к массиву <see cref="double"/></summary>
    /// <param name="s">Сигнал</param>
    public static implicit operator double[](SamplesDigitalSignal s) => s._Samples;
}