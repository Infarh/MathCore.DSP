#nullable enable
using System.Numerics;

using MathCore.DI;
using MathCore.DSP.Filters;

using Microsoft.Extensions.DependencyInjection;

namespace WpfTest.Services;

/// <summary>Результаты анализа частотных и временных характеристик БИХ-фильтра</summary>
/// <param name="Frequencies">Частоты расчёта, Гц</param>
/// <param name="MagnitudeDb">АЧХ в децибелах</param>
/// <param name="PhaseDegrees">Непрерывная ФЧХ в градусах</param>
/// <param name="GroupDelay">Групповая задержка в отсчётах</param>
/// <param name="ImpulseResponse">Отсчёты импульсной характеристики</param>
/// <param name="Poles">Полюса передаточной функции в z-плоскости</param>
/// <param name="Zeros">Нули передаточной функции в z-плоскости</param>
/// <param name="Order">Порядок фильтра</param>
/// <param name="SamplingFrequency">Частота дискретизации, Гц</param>
public sealed record FilterAnalysis(
    double[] Frequencies,
    double[] MagnitudeDb,
    double[] PhaseDegrees,
    double[] GroupDelay,
    double[] ImpulseResponse,
    Complex[] Poles,
    Complex[] Zeros,
    int Order,
    double SamplingFrequency);

/// <summary>Расчёт характеристик БИХ-фильтра по коэффициентам полиномов</summary>
[Service(ServiceLifetime.Singleton)]
public sealed class FilterAnalysisService
{
    /// <summary>Нижняя граница отображаемой АЧХ, дБ</summary>
    public const double MagnitudeFloorDb = -160;

    /// <summary>Рассчитать характеристики фильтра</summary>
    /// <param name="Filter">Исследуемый фильтр</param>
    /// <param name="SamplingFrequency">Частота дискретизации, Гц</param>
    /// <param name="Points">Число точек частотной характеристики</param>
    /// <returns>Результаты анализа</returns>
    public FilterAnalysis Calculate(IIR Filter, double SamplingFrequency, int Points = 2048)
    {
        ArgumentNullException.ThrowIfNull(Filter);
        if (SamplingFrequency <= 0) throw new ArgumentOutOfRangeException(nameof(SamplingFrequency));
        if (Points < 16) throw new ArgumentOutOfRangeException(nameof(Points));

        var a = Filter.A.ToArray();
        var b = Filter.B.ToArray();
        var nyquist = SamplingFrequency / 2;

        var frequencies = new double[Points];
        var magnitude = new double[Points];
        var phase = new double[Points];

        for (var i = 0; i < Points; i++)
        {
            var f = nyquist * i / (Points - 1);
            frequencies[i] = f;
            var h = FrequencyResponse(a, b, f / SamplingFrequency);
            magnitude[i] = Math.Max(MagnitudeFloorDb, 20 * Math.Log10(h.Magnitude + 1e-15));
            phase[i] = h.Phase;
        }

        UnwrapPhase(phase);
        var group_delay = GetGroupDelay(phase, frequencies, magnitude, SamplingFrequency, Filter.Order);

        for (var i = 0; i < Points; i++)
            phase[i] = magnitude[i] < -80 ? double.NaN : phase[i] * 180 / Math.PI;

        // Разрыв линии на скачке фазы в нулях передачи
        for (var i = 1; i < Points; i++)
        {
            if (double.IsNaN(phase[i]) || double.IsNaN(phase[i - 1])) continue;
            if (Math.Abs(phase[i] - phase[i - 1]) <= 90) continue;

            phase[i] = double.NaN;
            if (i + 1 < Points) phase[i + 1] = double.NaN;
        }

        // На частоте Найквиста фаза действительного коэффициента передачи не определена
        if (Points > 1) phase[^1] = double.NaN;

        return new(
            frequencies,
            magnitude,
            phase,
            group_delay,
            GetImpulseResponse(Filter),
            GetRoots(a.Reverse().ToArray()),
            GetZeros(b.Reverse().ToArray(), a.Length - b.Length),
            Filter.Order,
            SamplingFrequency);
    }

    private static Complex[] GetZeros(double[] ReversedNumerator, int ExtraZeros)
    {
        var zeros = GetRoots(ReversedNumerator).ToList();
        for (var i = 0; i < ExtraZeros; i++) zeros.Add(Complex.Zero);
        return [.. zeros];
    }

    private static Complex FrequencyResponse(double[] A, double[] B, double NormalizedFrequency)
    {
        var z = Complex.Exp(new Complex(0, -2 * Math.PI * NormalizedFrequency));
        return Evaluate(B, z) / Evaluate(A, z);

        static Complex Evaluate(double[] Coefficients, Complex X)
        {
            var y = new Complex(Coefficients[^1], 0);
            for (var i = Coefficients.Length - 2; i >= 0; i--)
                y = y * X + Coefficients[i];
            return y;
        }
    }

    private static void UnwrapPhase(double[] Phase)    {
        const double two_pi = 2 * Math.PI;
        for (var i = 1; i < Phase.Length; i++)
        {
            var delta = Phase[i] - Phase[i - 1];
            if (delta > Math.PI) delta -= two_pi * Math.Ceiling((delta - Math.PI) / two_pi);
            else if (delta < -Math.PI) delta += two_pi * Math.Ceiling((-delta - Math.PI) / two_pi);
            Phase[i] = Phase[i - 1] + delta;
        }
    }

    private static double[] GetGroupDelay(double[] Phase, double[] Frequencies, double[] Magnitude, double SamplingFrequency, int Order)
    {
        var count = Phase.Length;
        var result = new double[count];
        var step = 2 * Math.PI * (Frequencies[1] - Frequencies[0]) / SamplingFrequency;
        var limit = 10d * (Order + 1);

        for (var i = 0; i < count; i++)
        {
            var left = i > 0 ? i - 1 : i;
            var right = i < count - 1 ? i + 1 : i;
            var delta_phase = Phase[right] - Phase[left];
            var delta_index = right - left;

            if (delta_index == 0 || Magnitude[i] < -80)
            {
                result[i] = double.NaN;
                continue;
            }

            var delay = -delta_phase / (delta_index * step);
            result[i] = double.IsFinite(delay) && Math.Abs(delay) <= limit ? delay : double.NaN;
        }

        return result;
    }

    private static double[] GetImpulseResponse(IIR Filter)
    {
        var count = Math.Clamp(8 * (Filter.Order + 1), 128, 1024);
        var state = new double[Filter.Order + 1];
        var impulse = new double[count];

        for (var i = 0; i < count; i++)
            impulse[i] = Filter.Process(i == 0 ? 1 : 0, state);

        return impulse;
    }

    /// <summary>Вычислить корни полинома методом Дюрана-Кернера</summary>
    /// <param name="Coefficients">Коэффициенты полинома от младшей степени к старшей</param>
    /// <returns>Массив комплексных корней</returns>
    private static Complex[] GetRoots(IReadOnlyList<double> Coefficients)
    {
        var n = Coefficients.Count - 1;
        while (n > 0 && Math.Abs(Coefficients[n]) < 1e-15) n--;
        if (n <= 0) return [];

        var a = new double[n + 1];
        var lead = Coefficients[n];
        for (var i = 0; i <= n; i++) a[i] = Coefficients[i] / lead;

        var radius = 0d;
        for (var i = 0; i < n; i++) radius = Math.Max(radius, Math.Abs(a[i]));
        radius = 1 + radius;

        var roots = new Complex[n];
        for (var i = 0; i < n; i++)
            roots[i] = Complex.FromPolarCoordinates(0.7 * radius, 2 * Math.PI * (i + 0.5) / n);

        const int max_iterations = 512;
        const double eps = 1e-13;

        for (var iteration = 0; iteration < max_iterations; iteration++)
        {
            var max_delta = 0d;

            for (var i = 0; i < n; i++)
            {
                var p = Evaluate(a, roots[i]);
                var d = Complex.One;
                for (var j = 0; j < n; j++)
                    if (i != j) d *= roots[i] - roots[j];

                if (d == Complex.Zero) d = new(1e-12, 1e-12);

                var next = roots[i] - p / d;
                max_delta = Math.Max(max_delta, Complex.Abs(next - roots[i]));
                roots[i] = next;
            }

            if (max_delta < eps) break;
        }

        return roots;

        static Complex Evaluate(double[] Coefficients, Complex X)
        {
            var y = new Complex(Coefficients[^1], 0);
            for (var i = Coefficients.Length - 2; i >= 0; i--)
                y = y * X + Coefficients[i];
            return y;
        }
    }
}
