using System.Collections.ObjectModel;

using static System.Array;

namespace MathCore.DSP.Filters;

/// <summary>Комплексный фильтр с конечной импульсной характеристикой для IQ-потока</summary>
public class ComplexFIR
{
    /// <summary>Импульсная характеристика</summary>
    private readonly double[] _ImpulseResponse;

    /// <summary>Состояние синфазного канала</summary>
    private readonly double[] _StateI;

    /// <summary>Состояние квадратурного канала</summary>
    private readonly double[] _StateQ;

    /// <summary>Импульсная характеристика</summary>
    public ReadOnlyCollection<double> ImpulseResponse => AsReadOnly(_ImpulseResponse);

    /// <summary>Порядок фильтра</summary>
    public int Order => _ImpulseResponse.Length - 1;

    /// <summary>Инициализация нового комплексного FIR-фильтра</summary>
    /// <param name="ImpulseResponse">Отсчёты импульсной характеристики фильтра</param>
    public ComplexFIR(double[] ImpulseResponse)
    {
        ArgumentNullException.ThrowIfNull(ImpulseResponse);
        if (ImpulseResponse.Length == 0)
            throw new ArgumentException("Размер массива импульсной характеристики должен быть больше 0", nameof(ImpulseResponse));

        _ImpulseResponse = ImpulseResponse;
        _StateI = new double[ImpulseResponse.Length];
        _StateQ = new double[ImpulseResponse.Length];
    }

    /// <summary>Обработать один комплексный отсчёт</summary>
    /// <param name="Sample">Входной комплексный отсчёт</param>
    /// <returns>Выходной комплексный отсчёт</returns>
    public Complex Process(Complex Sample)
    {
        var i_filtered = _StateI.FilterSample(_ImpulseResponse, Sample.Re);
        var q_filtered = _StateQ.FilterSample(_ImpulseResponse, Sample.Im);

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
    public Complex FrequencyResponse(double f) => _ImpulseResponse.FrequencyResponse(f);
}
