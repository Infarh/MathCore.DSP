using System.Numerics;

using MathCore.DSP.Infrastructure;

namespace MathCore.DSP;

/// <summary>Методы-расширения для вещественных массивов</summary>
public static class DoubleArrayDSPExtensions
{
    internal static Vector<double> ToVector(this double[] array) => new(array);

    /// <summary>Вычислить значение коэффициента передачи фильтра, заданного импульсной характеристикой</summary>
    /// <param name="ImpulseResponse">Массив отсчётов импульсной характеристики</param>
    /// <param name="f">Частота вычисления коэффициента передачи</param>
    /// <param name="dt">Период дискретизации импульсной характеристики</param>
    /// <returns>Комплексное значение коэффициента передачи фильтра с указанной импульсной характеристикой</returns>
    public static Complex FrequencyResponse(this double[] ImpulseResponse, double f, double dt)
        => ImpulseResponse.FrequencyResponse(f * dt);

    /// <summary>Вычислить значение коэффициента передачи фильтра, заданного импульсной характеристикой</summary>
    /// <param name="ImpulseResponse">Массив отсчётов импульсной характеристики</param>
    /// <param name="f">Нормированная частота вычисления коэффициента передачи</param>
    /// <returns>Комплексное значение коэффициента передачи фильтра с указанной импульсной характеристикой</returns>
    public static Complex FrequencyResponse(this double[] ImpulseResponse, double f)
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
    public static double FilterSample(this double[] State, double[] ImpulseResponse, double Sample)
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
    public static double FilterSampleVector(this double[] State, double[] ImpulseResponse, double Sample)
    {
        Array.Copy(State, 0, State, 1, State.Length - 1);
        State[0] = Sample;

        return Vector.Dot(State.ToVector() * ImpulseResponse.ToVector(), Vector<double>.One);
    }

    public static IEnumerable<double> FilterFIR(
        this IEnumerable<double> samples,
        double[] ImpulseResponse,
        double[] State)
    {
        if (samples is null) throw new ArgumentNullException(nameof(samples));
        if (ImpulseResponse is null) throw new ArgumentNullException(nameof(ImpulseResponse));
        if (State is null) throw new ArgumentNullException(nameof(State));
        if (ImpulseResponse.Length != State.Length) throw new InvalidOperationException("Размер массива импульсной характеристики не соответствует размеру массива состояния фильтра");

        foreach (var sample in samples)
            yield return State.FilterSample(ImpulseResponse, sample);
    }

    public static IEnumerable<double> FilterFIR(this IEnumerable<double> samples, double[] ImpulseResponse)
        => samples.NotNull().FilterFIR(ImpulseResponse.NotNull(), new double[ImpulseResponse.Length]);


    public static Complex FrequencyResponse(double[] A, double[] B, double f, double dt)
        => FrequencyResponse(A.NotNull(), B.NotNull(), f * dt);

    public static Complex FrequencyResponse(this (IReadOnlyList<double> A, IReadOnlyList<double> B) Filter, double f, double dt) =>
        FrequencyResponse(Filter.A, Filter.B, f * dt);

    public static Complex FrequencyResponse(this (IReadOnlyList<double> A, IReadOnlyList<double> B) Filter, double f) =>
        FrequencyResponse(Filter.A, Filter.B, f);

    /// <summary>Расчёт коэффициента передачи рекуррентного фильтра, заданного массивами своих коэффициентов для указанной частоты</summary>
    /// <param name="A">Массив коэффициентов обратных связей</param>
    /// <param name="B">Массив коэффициентов прямых связей</param>
    /// <param name="f">Частота, на которой требуется рассчитать коэффициент передачи фильтра</param>
    /// <returns>Значение комплексного коэффициента передачи рекуррентного фильтра на заданной частоте</returns>
    public static Complex FrequencyResponse(IReadOnlyList<double> A, IReadOnlyList<double> B, double f)
    {
        var p = Complex.Exp(-Consts.pi2 * f);

        static Complex Sum(IReadOnlyList<double> V, Complex p)
        {
            var (re, im)     = (V[^1], 0d);
            var (e_re, e_im) = p;

            for (var i = V.Count - 2; i >= 0; i--) 
                (re, im) = (re * e_re - im * e_im + V[i], re * e_im + im * e_re);
            
            return new Complex(re, im);
        }

        return Sum(B, p) / Sum(A, p);
    }

    public static Complex DigitalFrequencyResponseFromZPoles(
        IEnumerable<Complex> ZerosZ,
        IEnumerable<Complex> PolesZ,
        double f,
        double dt)
    {
        var z = Complex.Exp(-Consts.pi2 * f * dt);

        var P0 = Complex.Real;
        var one = Complex.Real;
        foreach (var z0 in ZerosZ)
        {
            var zz = z0 * z;
            if (zz == one)
                return 0;
            P0 *= 1 - zz;
        }

        var Pp = Complex.Real;
        foreach (var zp in PolesZ)
        {
            var zz = zp * z;
            if (zz == one)
                return new Complex(double.PositiveInfinity, double.PositiveInfinity);
            Pp *= 1 - zz;
        }

        return P0 / Pp;
    }

    public static Complex AnalogFrequencyResponseFromPoles(
        IEnumerable<Complex> P0,
        IEnumerable<Complex> Pp,
        double f)
    {
        var p = Complex.ImValue(Consts.pi2 * f);

        var zeros = Complex.Real;
        foreach (var p0 in P0)
        {
            if (p0 == p)
                return 0;
            zeros *= p - p0;
        }

        var poles = Complex.Real;
        foreach (var pp in Pp)
        {
            if (p == pp)
                return new Complex(double.PositiveInfinity, double.PositiveInfinity);
            poles *= p - pp;
        }

        return zeros / poles;
    }

    /// <summary>Выполнение фильтрации очередного отсчёта цифрового сигнала с помощью коэффициентов рекуррентного фильтра</summary>
    /// <param name="State">Вектор состояния фильтра</param>
    /// <param name="A">Вектор коэффициентов обратных связей</param>
    /// <param name="B">Вектор коэффициентов прямых связей</param>
    /// <param name="Sample">Фильтруемый отсчёт</param>
    /// <returns>Обработанное значение</returns>
    public static double FilterSample(
        this double[] State,
        double[] A,
        double[] B,
        double Sample)
    {
        var a0 = 1 / A[0];

        var result = 0d;
        var input = Sample;
        var b_length = B.Length;
        if (A.Length == b_length)
            for (var i = State.Length - 1; i >= 1; i--)
            {
                //(State[i], result, input) = (State[i - 1], result + State[i - 1] * B[i] * a0, input - State[i - 1] * A[i] * a0);
                var v = State[i - 1];
                State[i] = v;
                result += v * B[i] * a0;
                input -= v * A[i] * a0;
                //(State[i], result, input) = (v, result + v * B[i] * a0, input - v * A[i] * a0);
            }
        else
        {
            for (var i = State.Length - 1; i >= b_length; i--)
            {
                var v = State[i - 1];
                State[i] = v;
                input -= v * A[i] * a0;
            }
            for (var i = b_length - 1; i >= 1; i--)
            {
                var v = State[i - 1];
                State[i] = v;
                result += v * B[i] * a0;
                input -= v * A[i] * a0;
            }
        }

        State[0] = input;
        return result + input * B[0] * a0;
    }

    public static IEnumerable<double> FilterIIR(
        this IEnumerable<double> samples,
        double[] A,
        double[] B,
        double[] State)
    {
        if (samples is null) throw new ArgumentNullException(nameof(samples));
        if (A is null) throw new ArgumentNullException(nameof(A));
        if (B is null) throw new ArgumentNullException(nameof(B));
        if (A.Length < B.Length) throw new InvalidOperationException("Размеры массивов числителя и знаменателя передаточной функции не равны");
        if (State is null) throw new ArgumentNullException(nameof(State));

        foreach (var sample in samples)
            yield return FilterSample(State, A, B, sample);
    }

    public static IEnumerable<double> FilterIIR(
        this IEnumerable<double> samples,
        double[] A,
        double[] B)
        => samples.FilterIIR(A, B, new double[A.Length]);
}