using static MathCore.Polynom.Array;

namespace MathCore.DSP.Filters;

public static partial class FilterEx
{
    /// <summary>Целочисленные коэффициенты IIR-фильтра в прямой форме для арифметики с фиксированной точкой</summary>
    public readonly record struct FixedPointIirCoefficients
    {
        /// <summary>Коэффициенты числителя</summary>
        public long[] B { get; }

        /// <summary>Коэффициенты знаменателя</summary>
        public long[] A { get; }

        /// <summary>Разрядность знакового целого коэффициента</summary>
        public int BitDepth { get; }

        /// <summary>Количество дробных битов</summary>
        public int FractionalBits { get; }

        /// <summary>Количество битов целой части без учёта знакового бита</summary>
        public int IntegerBits => BitDepth - FractionalBits - 1;

        /// <summary>Шаг квантования коэффициентов</summary>
        public double Scale => Math.Pow(2, -FractionalBits);

        /// <summary>Максимальная абсолютная ошибка квантования коэффициентов</summary>
        public double MaxError { get; }

        /// <summary>Инициализация структуры коэффициентов IIR-фильтра в прямой форме</summary>
        /// <param name="B">Коэффициенты числителя</param>
        /// <param name="A">Коэффициенты знаменателя</param>
        /// <param name="BitDepth">Разрядность коэффициентов</param>
        /// <param name="FractionalBits">Количество дробных битов</param>
        /// <param name="MaxError">Максимальная абсолютная ошибка квантования</param>
        public FixedPointIirCoefficients(long[] B, long[] A, int BitDepth, int FractionalBits, double MaxError)
        {
            this.B = B ?? throw new ArgumentNullException(nameof(B));
            this.A = A ?? throw new ArgumentNullException(nameof(A));
            this.BitDepth = BitDepth;
            this.FractionalBits = FractionalBits;
            this.MaxError = MaxError;
        }
    }

    /// <summary>Целочисленные коэффициенты одного SOS-звена для арифметики с фиксированной точкой</summary>
    public readonly record struct FixedPointSosSection
    {
        /// <summary>Коэффициенты числителя звена</summary>
        public long[] B { get; }

        /// <summary>Коэффициенты знаменателя звена</summary>
        public long[] A { get; }

        /// <summary>Разрядность знакового целого коэффициента</summary>
        public int BitDepth { get; }

        /// <summary>Количество дробных битов звена</summary>
        public int FractionalBits { get; }

        /// <summary>Количество битов целой части без учёта знакового бита</summary>
        public int IntegerBits => BitDepth - FractionalBits - 1;

        /// <summary>Шаг квантования коэффициентов звена</summary>
        public double Scale => Math.Pow(2, -FractionalBits);

        /// <summary>Максимальная абсолютная ошибка квантования коэффициентов звена</summary>
        public double MaxError { get; }

        /// <summary>Инициализация структуры коэффициентов SOS-звена</summary>
        /// <param name="B">Коэффициенты числителя звена</param>
        /// <param name="A">Коэффициенты знаменателя звена</param>
        /// <param name="BitDepth">Разрядность коэффициентов</param>
        /// <param name="FractionalBits">Количество дробных битов</param>
        /// <param name="MaxError">Максимальная абсолютная ошибка квантования</param>
        public FixedPointSosSection(long[] B, long[] A, int BitDepth, int FractionalBits, double MaxError)
        {
            this.B = B ?? throw new ArgumentNullException(nameof(B));
            this.A = A ?? throw new ArgumentNullException(nameof(A));
            this.BitDepth = BitDepth;
            this.FractionalBits = FractionalBits;
            this.MaxError = MaxError;
        }
    }

    /// <summary>Целочисленные коэффициенты IIR-фильтра в виде каскада SOS-звеньев</summary>
    public readonly record struct FixedPointIirSosCoefficients
    {
        /// <summary>Массив квантованных SOS-звеньев</summary>
        public FixedPointSosSection[] Sections { get; }

        /// <summary>Разрядность знакового целого коэффициента</summary>
        public int BitDepth { get; }

        /// <summary>Максимальная абсолютная ошибка квантования среди всех звеньев</summary>
        public double MaxError { get; }

        /// <summary>Инициализация структуры SOS-представления IIR-фильтра</summary>
        /// <param name="Sections">Массив квантованных SOS-звеньев</param>
        /// <param name="BitDepth">Разрядность коэффициентов</param>
        /// <param name="MaxError">Максимальная абсолютная ошибка квантования</param>
        public FixedPointIirSosCoefficients(FixedPointSosSection[] Sections, int BitDepth, double MaxError)
        {
            this.Sections = Sections ?? throw new ArgumentNullException(nameof(Sections));
            this.BitDepth = BitDepth;
            this.MaxError = MaxError;
        }
    }

    /// <summary>Экспортировать коэффициенты фильтра в целочисленную прямую форму для арифметики с фиксированной точкой</summary>
    /// <param name="Filter">Экспортируемый фильтр</param>
    /// <param name="BitDepth">Разрядность знакового целого коэффициента</param>
    /// <param name="FractionalBits">Число дробных битов или <see langword="null"/> для автоматического подбора</param>
    /// <returns>Квантованные коэффициенты прямой формы</returns>
    /// <exception cref="ArgumentNullException">Если <paramref name="Filter"/> равен <see langword="null"/></exception>
    /// <exception cref="NotSupportedException">Если фильтр не является <see cref="IIR"/></exception>
    public static FixedPointIirCoefficients ExportToFixedPoint(this Filter Filter, int BitDepth = 16, int? FractionalBits = null)
    {
        if (Filter is null) throw new ArgumentNullException(nameof(Filter));

        return Filter is IIR iir
            ? iir.ExportToFixedPoint(BitDepth, FractionalBits)
            : throw new NotSupportedException($"Экспорт в фиксированную точку поддерживается только для фильтров типа {nameof(IIR)}");
    }

    /// <summary>Экспортировать коэффициенты IIR-фильтра в целочисленную прямую форму для арифметики с фиксированной точкой</summary>
    /// <param name="Filter">Экспортируемый IIR-фильтр</param>
    /// <param name="BitDepth">Разрядность знакового целого коэффициента</param>
    /// <param name="FractionalBits">Число дробных битов или <see langword="null"/> для автоматического подбора</param>
    /// <returns>Квантованные коэффициенты прямой формы</returns>
    /// <exception cref="ArgumentNullException">Если <paramref name="Filter"/> равен <see langword="null"/></exception>
    public static FixedPointIirCoefficients ExportToFixedPoint(this IIR Filter, int BitDepth = 16, int? FractionalBits = null)
    {
        if (Filter is null) throw new ArgumentNullException(nameof(Filter));

        ValidateBitDepth(BitDepth);

        var (b_normalized, a_normalized) = Normalize(Filter.B, Filter.A);
        var coefficients = a_normalized.Concat(b_normalized).ToArray();
        var fractional_bits = ResolveFractionalBits(coefficients, BitDepth, FractionalBits);
        var a = Quantize(a_normalized, BitDepth, fractional_bits);
        var b = Quantize(b_normalized, BitDepth, fractional_bits);
        var max_error = Math.Max(
            GetMaxError(a_normalized, a, fractional_bits),
            GetMaxError(b_normalized, b, fractional_bits));

        return new(b, a, BitDepth, fractional_bits, max_error);
    }

    /// <summary>Экспортировать коэффициенты фильтра в целочисленное SOS-представление для арифметики с фиксированной точкой</summary>
    /// <param name="Filter">Экспортируемый фильтр</param>
    /// <param name="BitDepth">Разрядность знакового целого коэффициента</param>
    /// <param name="FractionalBits">Число дробных битов или <see langword="null"/> для автоматического подбора по каждому звену</param>
    /// <returns>Квантованные коэффициенты каскада SOS-звеньев</returns>
    /// <exception cref="ArgumentNullException">Если <paramref name="Filter"/> равен <see langword="null"/></exception>
    /// <exception cref="NotSupportedException">Если фильтр не является <see cref="IIR"/></exception>
    public static FixedPointIirSosCoefficients ExportToFixedPointSos(this Filter Filter, int BitDepth = 16, int? FractionalBits = null)
    {
        if (Filter is null) throw new ArgumentNullException(nameof(Filter));

        return Filter is IIR iir
            ? iir.ExportToFixedPointSos(BitDepth, FractionalBits)
            : throw new NotSupportedException($"Экспорт в SOS-представление поддерживается только для фильтров типа {nameof(IIR)}");
    }

    /// <summary>Экспортировать коэффициенты IIR-фильтра в целочисленное SOS-представление для арифметики с фиксированной точкой</summary>
    /// <param name="Filter">Экспортируемый IIR-фильтр</param>
    /// <param name="BitDepth">Разрядность знакового целого коэффициента</param>
    /// <param name="FractionalBits">Число дробных битов или <see langword="null"/> для автоматического подбора по каждому звену</param>
    /// <returns>Квантованные коэффициенты каскада SOS-звеньев</returns>
    /// <exception cref="ArgumentNullException">Если <paramref name="Filter"/> равен <see langword="null"/></exception>
    /// <exception cref="InvalidOperationException">Если не удалось выполнить секционное разложение фильтра</exception>
    public static FixedPointIirSosCoefficients ExportToFixedPointSos(this IIR Filter, int BitDepth = 16, int? FractionalBits = null)
    {
        if (Filter is null) throw new ArgumentNullException(nameof(Filter));

        ValidateBitDepth(BitDepth);

        var (b_normalized, a_normalized) = Normalize(Filter.B, Filter.A);
        var sections = BuildSos(a_normalized, b_normalized)
            .OrderBy(section => section.PoleRadius)
            .ToArray();
        var result_sections = new FixedPointSosSection[sections.Length];
        var max_error = 0d;

        for (var i = 0; i < sections.Length; i++)
        {
            var section = sections[i];
            var coefficients = section.A.Concat(section.B).ToArray();
            var fractional_bits = ResolveFractionalBits(coefficients, BitDepth, FractionalBits);
            var a = Quantize(section.A, BitDepth, fractional_bits);
            var b = Quantize(section.B, BitDepth, fractional_bits);
            var section_error = Math.Max(
                GetMaxError(section.A, a, fractional_bits),
                GetMaxError(section.B, b, fractional_bits));

            if (section_error > max_error)
                max_error = section_error;

            result_sections[i] = new(b, a, BitDepth, fractional_bits, section_error);
        }

        return new(result_sections, BitDepth, max_error);
    }

    private readonly record struct SosSection(double[] B, double[] A)
    {
        public double PoleRadius
        {
            get
            {
                var a1 = A[1];
                var a2 = A[2];
                var discriminant = Complex.Sqrt(a1 * a1 - 4 * a2);
                var root_1 = (-a1 + discriminant) / 2;
                var root_2 = (-a1 - discriminant) / 2;
                return Math.Max(root_1.Abs, root_2.Abs);
            }
        }
    }

    private static (double[] B, double[] A) Normalize(IReadOnlyList<double> B, IReadOnlyList<double> A)
    {
        if (B is null) throw new ArgumentNullException(nameof(B));
        if (A is null) throw new ArgumentNullException(nameof(A));
        if (B.Count == 0) throw new ArgumentException("Размер массива коэффициентов числителя должен быть больше 0", nameof(B));
        if (A.Count < 2) throw new ArgumentException("Размер массива коэффициентов знаменателя должен быть больше 1", nameof(A));

        var a0 = A[0];
        if (!IsFinite(a0) || Math.Abs(a0) <= double.Epsilon)
            throw new InvalidOperationException("Коэффициент A0 должен быть конечным и отличным от нуля");

        var normalized_b = B.Select(value => value / a0).ToArray();
        var normalized_a = A.Select(value => value / a0).ToArray();

        EnsureFinite(normalized_b, nameof(B));
        EnsureFinite(normalized_a, nameof(A));

        return (normalized_b, normalized_a);
    }

    private static void EnsureFinite(IEnumerable<double> Values, string ParamName)
    {
        if (Values is null) throw new ArgumentNullException(nameof(Values));
        if (ParamName is null) throw new ArgumentNullException(nameof(ParamName));

        foreach (var value in Values)
            if (!IsFinite(value))
                throw new InvalidOperationException($"Коэффициенты {ParamName} содержат неконечные значения");
    }

    private static bool IsFinite(double Value) => !double.IsNaN(Value) && !double.IsInfinity(Value);

    private static void ValidateBitDepth(int BitDepth)
    {
        if (BitDepth < 2 || BitDepth > 62)
            throw new ArgumentOutOfRangeException(nameof(BitDepth), BitDepth, "Разрядность коэффициентов должна находиться в диапазоне от 2 до 62 бит");
    }

    private static int ResolveFractionalBits(IReadOnlyList<double> Coefficients, int BitDepth, int? FractionalBits)
    {
        if (Coefficients is null) throw new ArgumentNullException(nameof(Coefficients));

        if (FractionalBits is { } explicit_fractional_bits)
        {
            ValidateFractionalBits(explicit_fractional_bits, BitDepth);
            ValidateQuantizationRange(Coefficients, BitDepth, explicit_fractional_bits);
            return explicit_fractional_bits;
        }

        var max_abs = Coefficients.Count == 0 ? 0d : Coefficients.Max(value => Math.Abs(value));
        if (!IsFinite(max_abs))
            throw new InvalidOperationException("Коэффициенты фильтра содержат некорректные значения");

        if (max_abs <= double.Epsilon)
            return BitDepth - 2;

        var max_value = GetMaxValue(BitDepth);
        var fractional_bits = (int)Math.Floor(Math.Log(max_value / max_abs, 2));
        if (fractional_bits < 0)
            fractional_bits = 0;
        if (fractional_bits > BitDepth - 2)
            fractional_bits = BitDepth - 2;

        while (fractional_bits >= 0)
        {
            if (FitsInRange(Coefficients, BitDepth, fractional_bits))
                return fractional_bits;

            fractional_bits--;
        }

        throw new InvalidOperationException("Коэффициенты фильтра не помещаются в заданную разрядность");
    }

    private static void ValidateFractionalBits(int FractionalBits, int BitDepth)
    {
        if (FractionalBits < 0 || FractionalBits > BitDepth - 2)
            throw new ArgumentOutOfRangeException(nameof(FractionalBits), FractionalBits, $"Число дробных битов должно находиться в диапазоне от 0 до {BitDepth - 2}");
    }

    private static void ValidateQuantizationRange(IReadOnlyList<double> Coefficients, int BitDepth, int FractionalBits)
    {
        if (!FitsInRange(Coefficients, BitDepth, FractionalBits))
            throw new InvalidOperationException("Коэффициенты фильтра не помещаются в заданную разрядность при указанном числе дробных битов");
    }

    private static bool FitsInRange(IReadOnlyList<double> Coefficients, int BitDepth, int FractionalBits)
    {
        var min_value = GetMinValue(BitDepth);
        var max_value = GetMaxValue(BitDepth);
        var scale = Math.Pow(2, FractionalBits);

        for (var i = 0; i < Coefficients.Count; i++)
        {
            var value = Coefficients[i];
            if (!IsFinite(value))
                return false;

            var quantized = Math.Round(value * scale, MidpointRounding.AwayFromZero);
            if (quantized < min_value || quantized > max_value)
                return false;
        }

        return true;
    }

    private static long[] Quantize(IReadOnlyList<double> Coefficients, int BitDepth, int FractionalBits)
    {
        var result = new long[Coefficients.Count];
        var scale = Math.Pow(2, FractionalBits);
        var min_value = GetMinValue(BitDepth);
        var max_value = GetMaxValue(BitDepth);

        for (var i = 0; i < Coefficients.Count; i++)
        {
            var value = Coefficients[i];
            var quantized = Math.Round(value * scale, MidpointRounding.AwayFromZero);
            if (quantized < min_value || quantized > max_value)
                throw new InvalidOperationException("Коэффициенты фильтра не помещаются в заданную разрядность");

            result[i] = (long)quantized;
        }

        return result;
    }

    private static double GetMaxError(IReadOnlyList<double> Source, IReadOnlyList<long> Quantized, int FractionalBits)
    {
        var scale = Math.Pow(2, -FractionalBits);
        var max_error = 0d;

        for (var i = 0; i < Source.Count; i++)
        {
            var error = Math.Abs(Source[i] - Quantized[i] * scale);
            if (error > max_error)
                max_error = error;
        }

        return max_error;
    }

    private static long GetMinValue(int BitDepth) => -(1L << (BitDepth - 1));

    private static long GetMaxValue(int BitDepth) => (1L << (BitDepth - 1)) - 1;

    private static SosSection[] BuildSos(double[] A, double[] B)
    {
        try
        {
            var poles = GetRoots(A);
            var zeros = GetRoots(B);
            var denominator_sections = BuildSectionsFromRoots(poles);
            if (denominator_sections.Length == 0)
                throw new InvalidOperationException("Не удалось сформировать секции знаменателя");

            var numerator_sections = BuildSectionsFromRoots(zeros);
            if (numerator_sections.Length > denominator_sections.Length)
                throw new InvalidOperationException("Число секций числителя не должно превышать число секций знаменателя");

            var sections = new SosSection[denominator_sections.Length];
            for (var i = 0; i < denominator_sections.Length; i++)
            {
                var denominator = denominator_sections[i];
                var numerator = i < numerator_sections.Length ? numerator_sections[i] : new[] { 1d, 0d, 0d };
                sections[i] = new(numerator, denominator);
            }

            var gain = B[0];
            DistributeGain(sections, ref gain);
            return sections;
        }
        catch (Exception error)
        {
            throw new InvalidOperationException("Не удалось выполнить SOS-декомпозицию фильтра", error);
        }
    }

    private static Complex[] GetRoots(double[] Coefficients)
    {
        if (Coefficients is null) throw new ArgumentNullException(nameof(Coefficients));

        var n = Coefficients.Length - 1;
        while (n > 0 && Math.Abs(Coefficients[n]) < 1e-18)
            n--;

        if (n <= 0)
            return [];

        var coefficients = new double[n + 1];
        Array.Copy(Coefficients, coefficients, n + 1);

        var lead = coefficients[n];
        for (var i = 0; i <= n; i++)
            coefficients[i] /= lead;

        var roots = new Complex[n];
        var radius = 0.5;
        for (var i = 0; i < n; i++)
            roots[i] = Complex.Exp(radius, 2 * Math.PI * i / n);

        const int max_iterations = 256;
        const double eps = 1e-12;

        for (var iteration = 0; iteration < max_iterations; iteration++)
        {
            var max_delta = 0d;

            for (var i = 0; i < n; i++)
            {
                var polynomial = Evaluate(coefficients, roots[i]);
                var derivative = Complex.Real;
                for (var j = 0; j < n; j++)
                {
                    if (i == j) continue;
                    derivative *= roots[i] - roots[j];
                }

                if (derivative == Complex.Zero) continue;

                var next = roots[i] - polynomial / derivative;
                var delta = (next - roots[i]).Abs;
                if (delta > max_delta)
                    max_delta = delta;
                roots[i] = next;
            }

            if (max_delta < eps)
                break;
        }

        return roots;

        static Complex Evaluate(double[] Coefficients, Complex X)
        {
            Complex result = Coefficients[^1];
            for (var i = Coefficients.Length - 2; i >= 0; i--)
                result = result * X + Coefficients[i];
            return result;
        }
    }

    private static double[][] BuildSectionsFromRoots(Complex[] Roots)
    {
        if (Roots is null) throw new ArgumentNullException(nameof(Roots));
        if (Roots.Length == 0) return [];

        const double imag_eps = 1e-9;

        var complex_roots = new List<Complex>(Roots.Length);
        var real_roots = new List<double>(Roots.Length);

        foreach (var root in Roots)
            if (Math.Abs(root.Im) <= imag_eps)
                real_roots.Add(root.Re);
            else
                complex_roots.Add(root);

        var sections = new List<double[]>();

        while (complex_roots.Count > 0)
        {
            var root_1 = complex_roots[0];
            complex_roots.RemoveAt(0);

            var pair_index = -1;
            var min_delta = double.PositiveInfinity;
            var conjugate = root_1.ComplexConjugate;
            for (var i = 0; i < complex_roots.Count; i++)
            {
                var delta = (complex_roots[i] - conjugate).Abs;
                if (delta >= min_delta) continue;
                min_delta = delta;
                pair_index = i;
            }

            if (pair_index < 0)
            {
                real_roots.Add(root_1.Re);
                continue;
            }

            var root_2 = complex_roots[pair_index];
            complex_roots.RemoveAt(pair_index);

            var sum = root_1 + root_2;
            var multiply = root_1 * root_2;
            sections.Add([1, -sum.Re, multiply.Re]);
        }

        real_roots.Sort((x, y) => Math.Abs(y).CompareTo(Math.Abs(x)));
        for (var i = 0; i < real_roots.Count; i += 2)
        {
            var root_1 = real_roots[i];
            var root_2 = i + 1 < real_roots.Count ? real_roots[i + 1] : 0d;
            sections.Add([1, -(root_1 + root_2), root_1 * root_2]);
        }

        return sections.ToArray();
    }

    private static void DistributeGain(SosSection[] Sections, ref double Gain)
    {
        if (Sections is null) throw new ArgumentNullException(nameof(Sections));
        if (Sections.Length == 0 || Gain == 0 || double.IsNaN(Gain) || double.IsInfinity(Gain))
            return;

        var section_gain = Math.Pow(Math.Abs(Gain), 1d / Sections.Length);
        if (double.IsNaN(section_gain) || double.IsInfinity(section_gain) || section_gain == 0)
            return;

        var sign = Math.Sign(Gain);
        for (var i = 0; i < Sections.Length; i++)
        {
            var section = Sections[i];
            var scale = i == 0 ? section_gain * sign : section_gain;
            var numerator = section.B.ToArray();
            for (var k = 0; k < numerator.Length; k++)
                numerator[k] *= scale;

            Sections[i] = new(numerator, section.A);
        }

        Gain = 1;
    }
}
