#if NET8_0_OR_GREATER

using System.Numerics;

namespace MathCore.DSP;

public static partial class DoubleArrayDSPExtensions
{
    /// <summary>Вычислить значение коэффициента передачи фильтра, заданного импульсной характеристикой</summary>
    /// <param name="ImpulseResponse">Массив отсчётов импульсной характеристики</param>
    /// <param name="f">Нормированная частота вычисления коэффициента передачи</param>
    /// <returns>Комплексное значение коэффициента передачи фильтра с указанной импульсной характеристикой</returns>
    public static Complex FrequencyResponse(Span<double> ImpulseResponse, double f)
    {
        var e = Complex.Exp(-Consts.pi2 * f);
        Complex result = ImpulseResponse[^1];
        for (var i = ImpulseResponse.Length - 2; i >= 0; i--)
            result = result * e + ImpulseResponse[i];
        return result;
    }

    /// <summary>Вычислить значение коэффициента передачи фильтра, заданного импульсной характеристикой</summary>
    /// <param name="ImpulseResponse">Массив отсчётов импульсной характеристики</param>
    /// <param name="f">Нормированная частота вычисления коэффициента передачи</param>
    /// <returns>Комплексное значение коэффициента передачи фильтра с указанной импульсной характеристикой</returns>
    public static Complex FrequencyResponse(ReadOnlySpan<double> ImpulseResponse, double f)
    {
        var e = Complex.Exp(-Consts.pi2 * f);
        Complex result = ImpulseResponse[^1];
        for (var i = ImpulseResponse.Length - 2; i >= 0; i--)
            result = result * e + ImpulseResponse[i];
        return result;
    }

    /// <summary>Вычисление выходного значения фильтра, заданного вектором состояния и импульсной характеристикой</summary>
    /// <param name="State">Вектор состояния фильтра</param>
    /// <param name="ImpulseResponse">Массив значений импульсной характеристики</param>
    /// <param name="Sample">Значение входного отсчёта фильтра</param>
    /// <returns>Значение выходного отсчёта фильтра</returns>
    public static double FilterSample(Span<double> State, Span<double> ImpulseResponse, double Sample)
    {
        var result = 0d;

        for (var i = State.Length - 1; i >= 1; i--)
        {
            State[i] = State[i - 1];
            result += State[i] * ImpulseResponse[i];
        }

        State[0] = Sample;

        return result + Sample * ImpulseResponse[0];
    }

    /// <summary>Вычисление выходного значения фильтра, заданного вектором состояния и импульсной характеристикой</summary>
    /// <param name="State">Вектор состояния фильтра</param>
    /// <param name="ImpulseResponse">Массив значений импульсной характеристики</param>
    /// <param name="Sample">Значение входного отсчёта фильтра</param>
    /// <returns>Значение выходного отсчёта фильтра</returns>
    public static double FilterSample(Span<double> State, ReadOnlySpan<double> ImpulseResponse, double Sample)
    {
        var result = 0d;

        for (var i = State.Length - 1; i >= 1; i--)
        {
            State[i] = State[i - 1];
            result += State[i] * ImpulseResponse[i];
        }

        State[0] = Sample;

        return result + Sample * ImpulseResponse[0];
    }

    /// <summary>Вычисление выходного значения фильтра, заданного вектором состояния и импульсной характеристикой</summary>
    /// <param name="State">Вектор состояния фильтра</param>
    /// <param name="ImpulseResponse">Массив значений импульсной характеристики</param>
    /// <param name="Sample">Значение входного отсчёта фильтра</param>
    /// <returns>Значение выходного отсчёта фильтра</returns>
    public static double FilterSampleVector(Span<double> State, Span<double> ImpulseResponse, double Sample)
    {
        State[..^1].CopyTo(State[1..]);
        State[0] = Sample;

        return Vector.Dot(new Vector<double>(State) * new Vector<double>(ImpulseResponse), Vector<double>.One);
    }
}

#endif