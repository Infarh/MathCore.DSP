using System.Collections.ObjectModel;

using static System.Array;

namespace MathCore.DSP.Filters;

/// <summary>Комплексный фильтр с бесконечной импульсной характеристикой для IQ-потока</summary>
public class ComplexIIR
{
    /// <summary>Коэффициенты прямой связи</summary>
    private readonly double[] _B;

    /// <summary>Коэффициенты обратной связи</summary>
    private readonly double[] _A;

    /// <summary>Состояние синфазного канала</summary>
    private readonly double[] _StateI;

    /// <summary>Состояние квадратурного канала</summary>
    private readonly double[] _StateQ;

    /// <summary>Массив коэффициентов полинома числителя</summary>
    public ReadOnlyCollection<double> B => AsReadOnly(_B);

    /// <summary>Массив коэффициентов полинома знаменателя</summary>
    public ReadOnlyCollection<double> A => AsReadOnly(_A);

    /// <summary>Порядок фильтра</summary>
    public int Order => _A.Length - 1;

    /// <summary>Инициализация нового комплексного IIR-фильтра</summary>
    /// <param name="B">Массив коэффициентов полинома числителя</param>
    /// <param name="A">Массив коэффициентов полинома знаменателя</param>
    public ComplexIIR(double[] B, double[] A)
    {
        ArgumentNullException.ThrowIfNull(B);
        ArgumentNullException.ThrowIfNull(A);

        if (B.Length == 0)
            throw new ArgumentException("Размер массива коэффициентов числителя должен быть больше 0", nameof(B));

        if (A.Length < 2)
            throw new ArgumentException("Размер массива коэффициентов знаменателя должен быть больше 1", nameof(A));

        if (B.Length > A.Length)
            throw new ArgumentException("Размер массива коэффициентов полинома числителя должен быть меньше, либо равен размеру массива коэффициентов полинома знаменателя");

        _B = B;
        _A = A;
        _StateI = new double[A.Length];
        _StateQ = new double[A.Length];
    }

    /// <summary>Обработать один комплексный отсчёт</summary>
    /// <param name="Sample">Входной комплексный отсчёт</param>
    /// <returns>Выходной комплексный отсчёт</returns>
    public Complex Process(Complex Sample)
    {
        var i_filtered = _StateI.FilterSample(_A, _B, Sample.Re);
        var q_filtered = _StateQ.FilterSample(_A, _B, Sample.Im);

        return new(i_filtered, q_filtered);
    }

    /// <summary>Обработать последовательность комплексных отсчётов</summary>
    /// <param name="Samples">Последовательность входных отсчётов</param>
    /// <returns>Последовательность выходных отсчётов</returns>
    public IEnumerable<Complex> Process(IEnumerable<Complex> Samples)
    {
        ArgumentNullException.ThrowIfNull(Samples);

        foreach (var sample in Samples)
            yield return Process(sample);
    }

    /// <summary>Обработать блок комплексных отсчётов</summary>
    /// <param name="Samples">Входной блок отсчётов</param>
    /// <param name="Destination">Выходной буфер</param>
    public void Process(ReadOnlySpan<Complex> Samples, Span<Complex> Destination)
    {
        if (Destination.Length < Samples.Length)
            throw new ArgumentException("Размер буфера назначения должен быть не меньше размера входного блока", nameof(Destination));

        for (var i = 0; i < Samples.Length; i++)
            Destination[i] = Process(Samples[i]);
    }

    /// <summary>Сбросить состояние фильтра</summary>
    public void Reset()
    {
        Clear(_StateI, 0, _StateI.Length);
        Clear(_StateQ, 0, _StateQ.Length);
    }

    /// <summary>Получить комплексный коэффициент передачи фильтра</summary>
    /// <param name="f">Нормированная частота</param>
    /// <returns>Значение коэффициента передачи</returns>
    public Complex FrequencyResponse(double f) => DoubleArrayDSPExtensions.FrequencyResponse(_A, _B, f);
}
