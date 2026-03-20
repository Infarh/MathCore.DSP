using System.Collections.ObjectModel;

using static System.Array;

using static MathCore.Polynom.Array;
// ReSharper disable ArgumentsStyleOther

namespace MathCore.DSP.Filters;

/// <summary>Фильтр с бесконечной импульсной характеристикой</summary>
public class IIR : DigitalFilter
{
    private readonly record struct SosSection(double B0, double B1, double B2, double A1, double A2);

    private const int SosOrderThreshold = 32;

    /// <summary>Массив коэффициентов полинома числителя</summary>
    private readonly double[] _B;
    /// <summary>Массив коэффициентов полинома знаменателя</summary>
    private readonly double[] _A;

    private readonly SosSection[]? _SosSections;
    private readonly double[]? _SosState;
    private readonly double _SosGain;

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

        if (Order >= SosOrderThreshold && TryCreateSos(A, B, out var sections, out var gain))
        {
            _SosSections = sections;
            _SosState = new double[sections.Length * 2];
            _SosGain = gain;
        }
    }

    public override double Process(double Sample, double[] state)
    {
        if (_SosSections is null || _SosState is null || !ReferenceEquals(state, State))
            return state.FilterSample(_A, _B, Sample);

        var x = Sample;
        for (var section_index = 0; section_index < _SosSections.Length; section_index++)
        {
            var section = _SosSections[section_index];
            var state_index = section_index * 2;

            var s1 = _SosState[state_index];
            var s2 = _SosState[state_index + 1];

            var y = section.B0 * x + s1;
            _SosState[state_index] = section.B1 * x - section.A1 * y + s2;
            _SosState[state_index + 1] = section.B2 * x - section.A2 * y;

            x = y;
        }

        return x * _SosGain;
    }

    public override double Process(double Sample) => Process(Sample, State);

    public override void Reset()
    {
        base.Reset();
        if (_SosState is not null)
            Clear(_SosState, 0, _SosState.Length);
    }

    public override Complex FrequencyResponse(double f) => DoubleArrayDSPExtensions.FrequencyResponse(_A, _B, f);

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

    private static bool TryCreateSos(double[] A, double[] B, out SosSection[] Sections, out double Gain)
    {
        Sections = [];
        Gain = 1;

        try
        {
            if (A.Length < 2 || Math.Abs(A[0]) <= double.Epsilon)
                return false;

            var a0 = A[0];
            var a_norm = A.ToArray(v => v / a0);
            var b_norm = B.ToArray(v => v / a0);

            var poles = GetRoots(a_norm);
            var zeros = GetRoots(b_norm);

            var den_sections = BuildSectionsFromRoots(poles);
            if (den_sections.Length == 0)
                return false;

            var num_sections = BuildSectionsFromRoots(zeros);
            if (num_sections.Length > den_sections.Length)
                return false;

            var sections = new SosSection[den_sections.Length];
            for (var i = 0; i < den_sections.Length; i++)
            {
                var (_, a1, a2) = den_sections[i];
                var (b0, b1, b2) = i < num_sections.Length ? num_sections[i] : (1d, 0d, 0d);
                sections[i] = new(b0, b1, b2, a1, a2);
            }

            var gain = b_norm[0];
            DistributeGain(sections, ref gain);

            if (!IsValidSos(sections, gain, a_norm, b_norm))
                return false;

            Sections = sections;
            Gain = gain;
            return true;
        }
        catch
        {
            return false;
        }
    }

    private static (double B0, double B1, double B2)[] BuildSectionsFromRoots(Complex[] Roots)
    {
        if (Roots.Length == 0)
            return [];

        const double imag_eps = 1e-9;

        var complex_roots = new List<Complex>(Roots.Length);
        var real_roots = new List<double>(Roots.Length);

        foreach (var root in Roots)
            if (Math.Abs(root.Im) <= imag_eps)
                real_roots.Add(root.Re);
            else
                complex_roots.Add(root);

        var sections = new List<(double B0, double B1, double B2)>();

        while (complex_roots.Count > 0)
        {
            var r1 = complex_roots[0];
            complex_roots.RemoveAt(0);

            var pair_index = -1;
            var min_delta = double.PositiveInfinity;
            var conjugate = r1.ComplexConjugate;
            for (var i = 0; i < complex_roots.Count; i++)
            {
                var delta = (complex_roots[i] - conjugate).Abs;
                if (delta >= min_delta) continue;
                min_delta = delta;
                pair_index = i;
            }

            if (pair_index < 0)
            {
                real_roots.Add(r1.Re);
                continue;
            }

            var r2 = complex_roots[pair_index];
            complex_roots.RemoveAt(pair_index);

            var sum = r1 + r2;
            var mult = r1 * r2;
            sections.Add((1, -sum.Re, mult.Re));
        }

        real_roots.Sort((x, y) => Math.Abs(y).CompareTo(Math.Abs(x)));
        for (var i = 0; i < real_roots.Count; i += 2)
        {
            var r1 = real_roots[i];
            var r2 = i + 1 < real_roots.Count ? real_roots[i + 1] : 0d;
            sections.Add((1, -(r1 + r2), r1 * r2));
        }

        return sections.ToArray();
    }

    private static void DistributeGain(SosSection[] Sections, ref double Gain)
    {
        if (Sections.Length == 0 || Gain == 0 || double.IsNaN(Gain) || double.IsInfinity(Gain))
            return;

        var section_gain = Math.Pow(Math.Abs(Gain), 1d / Sections.Length);
        if (double.IsNaN(section_gain) || double.IsInfinity(section_gain) || section_gain == 0)
            return;

        var sign = Math.Sign(Gain);
        for (var i = 0; i < Sections.Length; i++)
        {
            var s = Sections[i];
            var scale = i == 0 ? section_gain * sign : section_gain;
            Sections[i] = new(
                B0: s.B0 * scale,
                B1: s.B1 * scale,
                B2: s.B2 * scale,
                A1: s.A1,
                A2: s.A2);
        }

        Gain = 1;
    }

    private static bool IsValidSos(SosSection[] Sections, double Gain, double[] ANormalized, double[] BNormalized)
    {
        ArgumentNullException.ThrowIfNull(Sections);
        ArgumentNullException.ThrowIfNull(ANormalized);
        ArgumentNullException.ThrowIfNull(BNormalized);

        if (double.IsNaN(Gain) || double.IsInfinity(Gain))
            return false;

        const double max_abs = 1e6;
        const double max_pole_radius = 1.0001;

        foreach (var section in Sections)
        {
            if (double.IsNaN(section.B0) || double.IsInfinity(section.B0) ||
                double.IsNaN(section.B1) || double.IsInfinity(section.B1) ||
                double.IsNaN(section.B2) || double.IsInfinity(section.B2) ||
                double.IsNaN(section.A1) || double.IsInfinity(section.A1) ||
                double.IsNaN(section.A2) || double.IsInfinity(section.A2))
                return false;

            if (Math.Abs(section.B0) > max_abs || Math.Abs(section.B1) > max_abs || Math.Abs(section.B2) > max_abs ||
                Math.Abs(section.A1) > max_abs || Math.Abs(section.A2) > max_abs)
                return false;

            var discriminant = section.A1 * section.A1 - 4 * section.A2;
            if (discriminant >= 0)
            {
                var sqrt = Math.Sqrt(discriminant);
                var z1 = (-section.A1 + sqrt) / 2;
                var z2 = (-section.A1 - sqrt) / 2;
                if (Math.Abs(z1) > max_pole_radius || Math.Abs(z2) > max_pole_radius)
                    return false;
            }
            else
            {
                var re = -section.A1 / 2;
                var im = Math.Sqrt(-discriminant) / 2;
                var radius = Math.Sqrt(re * re + im * im);
                if (radius > max_pole_radius)
                    return false;
            }
        }

        var a_from_sections = new[] { 1d };
        var b_from_sections = new[] { Gain };
        foreach (var section in Sections)
        {
            a_from_sections = Multiply(a_from_sections, [1d, section.A1, section.A2]);
            b_from_sections = Multiply(b_from_sections, [section.B0, section.B1, section.B2]);
        }

        if (!IsPolynomialClose(a_from_sections, ANormalized))
            return false;

        if (!IsPolynomialClose(b_from_sections, BNormalized))
            return false;

        return true;
    }

    private static bool IsPolynomialClose(double[] Actual, double[] Expected)
    {
        const double abs_eps = 1e-8;
        const double rel_eps = 1e-4;

        if (Actual.Length != Expected.Length)
            return false;

        for (var i = 0; i < Actual.Length; i++)
        {
            var actual = Actual[i];
            var expected = Expected[i];
            var delta = Math.Abs(actual - expected);
            var scale = Math.Max(1, Math.Abs(expected));
            if (delta <= abs_eps) continue;
            if (delta / scale <= rel_eps) continue;
            return false;
        }

        return true;
    }

    private static Complex[] GetRoots(double[] Coefficients)
    {
        var n = Coefficients.Length - 1;
        while (n > 0 && Math.Abs(Coefficients[n]) < 1e-18)
            n--;

        if (n <= 0)
            return [];

        var c = new double[n + 1];
        Array.Copy(Coefficients, c, n + 1);

        var lead = c[n];
        for (var i = 0; i <= n; i++)
            c[i] /= lead;

        var roots = new Complex[n];
        var radius = 0.5;
        for (var i = 0; i < n; i++)
            roots[i] = Complex.Exp(radius, Consts.pi2 * i / n);

        const int max_iterations = 256;
        const double eps = 1e-12;

        for (var iteration = 0; iteration < max_iterations; iteration++)
        {
            var max_delta = 0d;

            for (var i = 0; i < n; i++)
            {
                var p = Evaluate(c, roots[i]);
                var d = Complex.Real;
                for (var j = 0; j < n; j++)
                {
                    if (i == j) continue;
                    d *= roots[i] - roots[j];
                }

                if (d == Complex.Zero) continue;

                var next = roots[i] - p / d;
                var delta = (next - roots[i]).Abs;
                if (delta > max_delta) max_delta = delta;
                roots[i] = next;
            }

            if (max_delta < eps)
                break;
        }

        return roots;

        static Complex Evaluate(double[] c, Complex x)
        {
            Complex y = c[^1];
            for (var i = c.Length - 2; i >= 0; i--)
                y = y * x + c[i];
            return y;
        }
    }
}