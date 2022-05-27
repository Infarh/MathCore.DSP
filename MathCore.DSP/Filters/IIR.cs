using System.Collections.ObjectModel;

using static System.Array;

using static MathCore.Polynom.Array;
// ReSharper disable ArgumentsStyleOther

namespace MathCore.DSP.Filters;

/// <summary>Фильтр с бесконечной импульсной характеристикой</summary>
public class IIR : DigitalFilter
{
    /// <summary>Массив коэффициентов полинома числителя</summary>
    private readonly double[] _B;
    /// <summary>Массив коэффициентов полинома знаменателя</summary>
    private readonly double[] _A;

    /// <summary>Массив коэффициентов полинома числителя</summary>
    public ReadOnlyCollection<double> B => AsReadOnly(_B);

    /// <summary>Массив коэффициентов полинома знаменателя</summary>
    public ReadOnlyCollection<double> A => AsReadOnly(_A);

    /// <summary>Инициализация нового цифрового фильтра с бесконечной импульсной характеристикой</summary>
    /// <param name="B">Массив коэффициентов полинома числителя</param>
    /// <param name="A">Массив коэффициентов полинома знаменателя</param>
    /// <exception cref="ArgumentException">Если число коэффициентов полинома числителя == 0</exception>
    /// <exception cref="ArgumentException">Если число коэффициентов полинома знаменателя меньше 2</exception>
    /// <exception cref="ArgumentException">Число коэффициентов полинома числителя должно быть меньше числа коэффициентов знаменателя</exception>
    public IIR(double[] B, double[] A)
        : base(
            Math.Max(
                (A ?? throw new ArgumentNullException(nameof(A))).Length,
                (B ?? throw new ArgumentNullException(nameof(B))).Length))
    {
        if (B.Length == 0) throw new ArgumentException("Размер массива коэффициентов числителя должен быть больше 0", nameof(B));
        if (A.Length < 2) throw new ArgumentException("Размер массива коэффициентов знаменателя должен быть больше 1", nameof(A));
        if (B.Length > A.Length) throw new ArgumentException("Размер массива коэффициентов полинома числителя должен быть меньше, либо равен размеру массива коэффициентов полинома знаменателя");
        _B = B;
        _A = A;
    }

    public override double Process(double Sample, double[] state) => state.FilterSample(_A, _B, Sample);

    public override double Process(double Sample) => Process(Sample, State);

    public override Complex GetTransmissionCoefficient(double f) => DoubleArrayDSPExtensions.FrequencyResponse(_A, _B, f);

    /// <summary>Последовательное соединение фильтра с другим <see cref="IIR"/></summary>
    /// <param name="filter">Соединяемый фильтр</param>
    /// <returns>Фильтр, представляющий собой результат последовательного соединения двух фильтров</returns>
    public IIR ConnectionSerialTo(IIR filter) =>
        filter is null
            ? throw new ArgumentNullException(nameof(filter))
            : new IIR(
                B: Multiply(_B, filter._B),
                A: Multiply(_A, filter._A));

    /// <summary>Параллельное соединение фильтра с другим <see cref="IIR"/></summary>
    /// <param name="filter">Соединяемый фильтр</param>
    /// <returns>Фильтр, представляющий собой результат параллельного соединения двух фильтров</returns>
    public IIR ConnectionParallelTo(IIR filter) =>
        filter is null
            ? throw new ArgumentNullException(nameof(filter))
            : new IIR(
                B: Multiply(
                    Sum(
                        Multiply(_B, filter._A), 
                        Multiply(filter._B, _A)), 
                    0.5),
                A: Multiply(_A, filter._A));
}