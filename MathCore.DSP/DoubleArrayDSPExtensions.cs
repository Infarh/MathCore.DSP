using System.Numerics;

namespace MathCore.DSP;

/// <summary>Методы-расширения для вещественных массивов</summary>
public static partial class DoubleArrayDSPExtensions
{
    /// <param name="array">Массив отсчётов импульсной характеристики</param>
    extension(double[] array)
    {
        internal Vector<double> ToVector() => new(array);

        /// <summary>Вычислить значение коэффициента передачи фильтра, заданного импульсной характеристикой</summary>
        /// <param name="f">Частота вычисления коэффициента передачи</param>
        /// <param name="dt">Период дискретизации импульсной характеристики</param>
        /// <returns>Комплексное значение коэффициента передачи фильтра с указанной импульсной характеристикой</returns>
        public Complex FrequencyResponse(double f, double dt)
            => array.FrequencyResponse(f * dt);

        /// <summary>Вычислить значение коэффициента передачи фильтра, заданного импульсной характеристикой</summary>
        /// <param name="f">Нормированная частота вычисления коэффициента передачи</param>
        /// <returns>Комплексное значение коэффициента передачи фильтра с указанной импульсной характеристикой</returns>
        public Complex FrequencyResponse(double f)
        {
            var e = Complex.Exp(-Consts.pi2 * f);
            Complex result = array[^1];
            for (var i = array.Length - 2; i >= 0; i--)
                result = result * e + array[i];
            return result;
        }

        /// <summary>Вычисление выходного значения фильтра, заданного вектором состояния и импульсной характеристикой</summary>
        /// <param name="ImpulseResponse">Массив значений импульсной характеристики</param>
        /// <param name="Sample">Значение входного отсчёта фильтра</param>
        /// <returns>Значение выходного отсчёта фильтра</returns>
        public double FilterSample(double[] ImpulseResponse, double Sample)
        {
            var result = 0d;

            for (var i = array.Length - 1; i >= 1; i--)
            {
                array[i] = array[i - 1];
                result += array[i] * ImpulseResponse[i];
            }

            array[0] = Sample;

            return result + Sample * ImpulseResponse[0];
        }

        /// <summary>Вычисление выходного значения фильтра, заданного вектором состояния и импульсной характеристикой</summary>
        /// <param name="ImpulseResponse">Массив значений импульсной характеристики</param>
        /// <param name="Sample">Значение входного отсчёта фильтра</param>
        /// <returns>Значение выходного отсчёта фильтра</returns>
        public double FilterSampleVector(double[] ImpulseResponse, double Sample)
        {
            Array.Copy(array, 0, array, 1, array.Length - 1);
            array[0] = Sample;

            return Vector.Dot(array.ToVector() * ImpulseResponse.ToVector(), Vector<double>.One);
        }

        /// <summary>Выполнение фильтрации очередного отсчёта цифрового сигнала с помощью коэффициентов рекуррентного фильтра</summary>
        /// <param name="A">Вектор коэффициентов обратных связей</param>
        /// <param name="B">Вектор коэффициентов прямых связей</param>
        /// <param name="Sample">Фильтруемый отсчёт</param>
        /// <returns>Обработанное значение</returns>
        public double FilterSample(double[] A, double[] B, double Sample)
        {
            var a0 = 1 / A[0];

            var result = 0d;
            var input = Sample;
            var b_length = B.Length;
            if (A.Length == b_length)
                for (var i = array.Length - 1; i >= 1; i--)
                {
                    //(State[i], result, input) = (State[i - 1], result + State[i - 1] * B[i] * a0, input - State[i - 1] * A[i] * a0);
                    var v = array[i - 1];
                    array[i] = v;
                    result += v * B[i] * a0;
                    input -= v * A[i] * a0;
                    //(State[i], result, input) = (v, result + v * B[i] * a0, input - v * A[i] * a0);
                }
            else
            {
                for (var i = array.Length - 1; i >= b_length; i--)
                {
                    var v = array[i - 1];
                    array[i] = v;
                    input -= v * A[i] * a0;
                }
                for (var i = b_length - 1; i >= 1; i--)
                {
                    var v = array[i - 1];
                    array[i] = v;
                    result += v * B[i] * a0;
                    input -= v * A[i] * a0;
                }
            }

            array[0] = input;
            return result + input * B[0] * a0;
        }
    }

    /// <param name="samples">Последовательность отсчётов цифрового сигнала для фильтрации</param>
    extension(IEnumerable<double> samples)
    {
        /// <summary>Выполнение фильтрации цифрового сигнала с помощью FIR-фильтра</summary>
        /// <param name="ImpulseResponse">Импульсная характеристика фильтра</param>
        /// <param name="State">Вектор состояния фильтра</param>
        /// <returns>Перечисление отсчётов на выходе фильтра</returns>
        /// <exception cref="ArgumentNullException">Передана пустая ссылка в одном из параметров</exception>
        /// <exception cref="InvalidOperationException">Размер массива импульсной характеристики не соответствует размеру массива состояния фильтра</exception>
        public IEnumerable<double> FilterFIR(double[] ImpulseResponse, double[] State)
        {
            ArgumentNullException.ThrowIfNull(samples);
            ArgumentNullException.ThrowIfNull(ImpulseResponse);
            ArgumentNullException.ThrowIfNull(State);
            if (ImpulseResponse.Length != State.Length) 
                throw new InvalidOperationException("Размер массива импульсной характеристики не соответствует размеру массива состояния фильтра");

            foreach (var sample in samples)
                yield return State.FilterSample(ImpulseResponse, sample);
        }

        public IEnumerable<double> FilterFIR(double[] ImpulseResponse) => samples.NotNull().FilterFIR(ImpulseResponse.NotNull(), new double[ImpulseResponse.Length]);
    }

    public static Complex FrequencyResponse(double[] A, double[] B, double f, double dt)
        => FrequencyResponse(A.NotNull(), B.NotNull(), f * dt);

    extension((IReadOnlyList<double> A, IReadOnlyList<double> B) Filter)
    {
        public Complex FrequencyResponse(double f, double dt) => FrequencyResponse(Filter.A, Filter.B, f * dt);

        public Complex FrequencyResponse(double f) => FrequencyResponse(Filter.A, Filter.B, f);
    }

    /// <summary>Расчёт коэффициента передачи рекуррентного фильтра, заданного массивами своих коэффициентов для указанной частоты</summary>
    /// <param name="A">Массив коэффициентов обратных связей</param>
    /// <param name="B">Массив коэффициентов прямых связей</param>
    /// <param name="f">Частота, на которой требуется рассчитать коэффициент передачи фильтра</param>
    /// <returns>Значение комплексного коэффициента передачи рекуррентного фильтра на заданной частоте</returns>
    public static Complex FrequencyResponse(IReadOnlyList<double> A, IReadOnlyList<double> B, double f)
    {
        var p = Complex.Exp(-Consts.pi2 * f);

        return Sum(B, p) / Sum(A, p);

        static Complex Sum(IReadOnlyList<double> V, Complex p)
        {
            var (re, im)     = (V[^1], 0d);
            var (e_re, e_im) = p;

            for (var i = V.Count - 2; i >= 0; i--) 
                (re, im) = (re * e_re - im * e_im + V[i], re * e_im + im * e_re);
            
            return new(re, im);
        }
    }

    public static Complex DigitalFrequencyResponseFromZPoles(
        IEnumerable<Complex> ZerosZ,
        IEnumerable<Complex> PolesZ,
        double f,
        double dt)
    {
        var z = Complex.Exp(-Consts.pi2 * f * dt);

        var p0 = Complex.Real;
        var one = Complex.Real;
        foreach (var z0 in ZerosZ)
        {
            var zz = z0 * z;
            if (zz == one)
                return 0;
            p0 *= 1 - zz;
        }

        var pp = Complex.Real;
        foreach (var zp in PolesZ)
        {
            var zz = zp * z;
            if (zz == one)
                return new(double.PositiveInfinity, double.PositiveInfinity);
            pp *= 1 - zz;
        }

        return p0 / pp;
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
                return new(double.PositiveInfinity, double.PositiveInfinity);
            poles *= p - pp;
        }

        return zeros / poles;
    }

    extension(IEnumerable<double> samples)
    {
        public IEnumerable<double> FilterIIR(double[] A, double[] B, double[] State)
        {
            ArgumentNullException.ThrowIfNull(samples);
            ArgumentNullException.ThrowIfNull(A);
            ArgumentNullException.ThrowIfNull(B);
            if (A.Length < B.Length) throw new InvalidOperationException("Размеры массивов числителя и знаменателя передаточной функции не равны");
            ArgumentNullException.ThrowIfNull(State);

            foreach (var sample in samples)
                yield return FilterSample(State, A, B, sample);
        }

        public IEnumerable<double> FilterIIR(double[] A, double[] B) => samples.FilterIIR(A, B, new double[A.Length]);
    }
}