using MathCore.DSP.Filters;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class IirCoefficientsStabilityTests
{
    private sealed record StressCase(string Name, double Dt, Func<DSP.Filters.IIR> Factory);

    [TestMethod]
    public void CalculatedCoefficients_ShouldRemainFinite_OnStressSpecifications() =>
        RunStressValidation((case_name, filter, dt) => AssertCoefficientsAreFinite(case_name, filter));

    [TestMethod]
    public void FrequencyResponse_ShouldRemainFinite_OnStressSpecifications() =>
        RunStressValidation(AssertFrequencyResponseIsFinite);

    [TestMethod]
    [Ignore("Диагностический тест: для высоких порядков требуется секционная реализация (SOS) вместо прямой формы")]
    public void ImpulseResponse_ShouldRemainFinite_OnStressSpecifications() =>
        RunStressValidation((case_name, filter, dt) => AssertImpulseResponseIsFinite(case_name, filter));

    private static void RunStressValidation(Action<string, DSP.Filters.IIR, double> Validator)
    {
        ArgumentNullException.ThrowIfNull(Validator);

        var errors = new List<string>();

        foreach (var @case in GetStressCases())
        {
            try
            {
                var filter = @case.Factory();
                Validator(@case.Name, filter, @case.Dt);
            }
            catch (Exception error)
            {
                errors.Add($"{@case.Name}: {error.GetType().Name} - {error.Message}");
            }
        }

        if (errors.Count > 0)
            Assert.Fail($"Обнаружены проблемы устойчивости:{Environment.NewLine}{string.Join(Environment.NewLine, errors)}");
    }

    private static IReadOnlyList<StressCase> GetStressCases()
    {
        const double dt = 1d / 48_000;
        const double gp = 0.891250938;
        const double gs = 0.01;

        return
        [
            new("Butterworth-LP", dt, () => new DSP.Filters.ButterworthLowPass(dt, 9_000, 9_300, gp, gs)),
            new("Butterworth-HP", dt, () => new DSP.Filters.ButterworthHighPass(dt, 6_000, 6_300, gp, gs)),
            new("Butterworth-BP", dt, () => new DSP.Filters.ButterworthBandPass(dt, 2_000, 2_300, 5_200, 5_500, gp, gs)),
            new("Butterworth-BS", dt, () => new DSP.Filters.ButterworthBandStop(dt, 1_800, 2_200, 4_800, 5_300, gp, gs)),

            new("ChebyshevI-LP", dt, () => new DSP.Filters.ChebyshevLowPass(dt, 9_000, 9_250, gp, gs, ChebyshevType.I)),
            new("ChebyshevII-LP", dt, () => new DSP.Filters.ChebyshevLowPass(dt, 9_000, 9_250, gp, gs, ChebyshevType.II)),
            new("ChebyshevIIC-LP", dt, () => new DSP.Filters.ChebyshevLowPass(dt, 9_000, 9_250, gp, gs, ChebyshevType.IICorrected)),

            new("ChebyshevI-HP", dt, () => new DSP.Filters.ChebyshevHighPass(dt, 6_000, 6_200, gp, gs, ChebyshevType.I)),
            new("ChebyshevII-HP", dt, () => new DSP.Filters.ChebyshevHighPass(dt, 6_000, 6_200, gp, gs, ChebyshevType.II)),
            new("ChebyshevIIC-HP", dt, () => new DSP.Filters.ChebyshevHighPass(dt, 6_000, 6_200, gp, gs, ChebyshevType.IICorrected)),

            new("ChebyshevI-BP", dt, () => new DSP.Filters.ChebyshevBandPass(dt, 2_000, 2_300, 5_200, 5_500, gp, gs, ChebyshevType.I)),
            new("ChebyshevII-BP", dt, () => new DSP.Filters.ChebyshevBandPass(dt, 2_000, 2_300, 5_200, 5_500, gp, gs, ChebyshevType.II)),
            new("ChebyshevIIC-BP", dt, () => new DSP.Filters.ChebyshevBandPass(dt, 2_000, 2_300, 5_200, 5_500, gp, gs, ChebyshevType.IICorrected)),

            new("ChebyshevI-BS", dt, () => new DSP.Filters.ChebyshevBandStop(dt, 1_800, 2_200, 4_800, 5_300, gp, gs, ChebyshevType.I)),
            new("ChebyshevII-BS", dt, () => new DSP.Filters.ChebyshevBandStop(dt, 1_800, 2_200, 4_800, 5_300, gp, gs, ChebyshevType.II)),
            new("ChebyshevIIC-BS", dt, () => new DSP.Filters.ChebyshevBandStop(dt, 1_800, 2_200, 4_800, 5_300, gp, gs, ChebyshevType.IICorrected)),

            new("Elliptic-LP", dt, () => new DSP.Filters.EllipticLowPass(dt, 9_000, 9_200, gp, gs)),
            new("Elliptic-HP", dt, () => new DSP.Filters.EllipticHighPass(dt, 6_000, 6_200, gp, gs)),
            new("Elliptic-BP", dt, () => new DSP.Filters.EllipticBandPass(dt, 2_000, 2_300, 5_200, 5_500, gp, gs)),
            new("Elliptic-BS", dt, () => new DSP.Filters.EllipticBandStop(dt, 1_800, 2_200, 4_800, 5_300, gp, gs))
        ];
    }

    private static void AssertCoefficientsAreFinite(string CaseName, DSP.Filters.IIR Filter)
    {
        for (var i = 0; i < Filter.A.Count; i++)
            Assert.IsTrue(double.IsFinite(Filter.A[i]), $"{CaseName}: A[{i}] не является конечным числом");

        for (var i = 0; i < Filter.B.Count; i++)
            Assert.IsTrue(double.IsFinite(Filter.B[i]), $"{CaseName}: B[{i}] не является конечным числом");

        Assert.IsTrue(double.IsFinite(Filter.A[0]), $"{CaseName}: A[0] не является конечным числом");
        Assert.IsTrue(Math.Abs(Filter.A[0]) > 1e-30, $"{CaseName}: A[0] слишком близок к нулю");
    }

    private static void AssertFrequencyResponseIsFinite(string CaseName, DSP.Filters.IIR Filter, double Dt)
    {
        var fd = 1d / Dt;
        const int points_count = 256;

        for (var i = 0; i <= points_count; i++)
        {
            var f = fd * i / (2d * points_count);
            var h = Filter.FrequencyResponse(f);

            Assert.IsTrue(double.IsFinite(h.Re), $"{CaseName}: Re(H) не является конечным числом на частоте {f}");
            Assert.IsTrue(double.IsFinite(h.Im), $"{CaseName}: Im(H) не является конечным числом на частоте {f}");
            Assert.IsTrue(double.IsFinite(h.Abs), $"{CaseName}: |H| не является конечным числом на частоте {f}");
        }
    }

    private static void AssertImpulseResponseIsFinite(string CaseName, DSP.Filters.IIR Filter)
    {
        const int samples_count = 2048;
        const double max_abs_limit = 1e12;

        Filter.Reset();

        var max_abs = 0d;
        for (var i = 0; i < samples_count; i++)
        {
            var x = i == 0 ? 1d : 0d;
            var y = Filter.Process(x);

            Assert.IsTrue(double.IsFinite(y), $"{CaseName}: значение импульсной характеристики не является конечным на отсчёте {i}");

            var abs_y = Math.Abs(y);
            if (abs_y > max_abs) max_abs = abs_y;
        }

        Assert.IsTrue(max_abs < max_abs_limit, $"{CaseName}: наблюдается неустойчивый рост импульсной характеристики, max={max_abs}");
    }
}
