using MathCore.DSP.Filters;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class FilterExExportTests
{
    private const double SosFrequencyResponseTolerance = 6e-3;

    [TestMethod]
    public void ExportToFixedPoint_ShouldNormalizeCoefficients()
    {
        Filter filter = new IIR(
            B: [0.25, -0.5, 0.25],
            A: [2, -1, 0.5]);

        var coefficients = filter.ExportToFixedPoint(16);

        Assert.AreEqual(16, coefficients.BitDepth);
        Assert.AreEqual(14, coefficients.FractionalBits);
        CollectionAssert.AreEqual(new long[] { 2_048, -4_096, 2_048 }, coefficients.B);
        CollectionAssert.AreEqual(new long[] { 16_384, -8_192, 4_096 }, coefficients.A);
        Assert.AreEqual(0, coefficients.MaxError, 1e-15);
    }

    [TestMethod]
    public void ExportToFixedPointSos_ShouldPreserveFrequencyResponse()
    {
        var section_1 = new IIR(
            B: [1, -0.2, 0.04],
            A: [1, -1.2, 0.36]);
        var section_2 = new IIR(
            B: [1, 0.1, 0.01],
            A: [1, -1.6, 0.81]);
        Filter filter = section_1.ConnectionSerialTo(section_2);

        var coefficients = filter.ExportToFixedPointSos(16);

        Assert.HasCount(2, coefficients.Sections);

        foreach (var section in coefficients.Sections)
        {
            Assert.HasCount(3, section.A);
            Assert.HasCount(3, section.B);
            Assert.IsGreaterThan(0, section.A[0]);
        }

        var iir = (IIR)filter;
        foreach (var frequency in new[] { 0d, 0.05, 0.125, 0.2, 0.35 })
        {
            var expected = iir.FrequencyResponse(frequency);
            var actual = GetFrequencyResponse(coefficients, frequency);
            var delta = expected - actual;
            Assert.IsLessThan(SosFrequencyResponseTolerance, delta.Abs, $"Отклонение АЧХ слишком велико на частоте {frequency}: {delta.Abs}");
        }
    }

    private static Complex GetFrequencyResponse(FilterEx.FixedPointIirSosCoefficients Coefficients, double Frequency)
    {
        var response = new Complex(1, 0);
        foreach (var section in Coefficients.Sections)
            response *= GetFrequencyResponse(section, Frequency);

        return response;
    }

    private static Complex GetFrequencyResponse(FilterEx.FixedPointSosSection Section, double Frequency)
    {
        var angle = -2 * Math.PI * Frequency;
        var z1 = new Complex(Math.Cos(angle), Math.Sin(angle));
        var z2 = z1 * z1;
        var scale = Section.Scale;
        var numerator =
            Section.B[0] * scale +
            Section.B[1] * scale * z1 +
            Section.B[2] * scale * z2;
        var denominator =
            Section.A[0] * scale +
            Section.A[1] * scale * z1 +
            Section.A[2] * scale * z2;

        return numerator / denominator;
    }
}
