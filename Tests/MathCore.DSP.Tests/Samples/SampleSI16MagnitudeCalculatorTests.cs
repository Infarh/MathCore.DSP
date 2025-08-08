using MathCore.DSP.Samples;

namespace MathCore.DSP.Tests.Samples;

[TestClass]
public class SampleSI16MagnitudeCalculatorTests
{
    /// <summary>Проверяет корректность вычисления модуля для всех комбинаций I и Q типа sbyte</summary>
    [TestMethod]
    public void GetMagnitude_AllSByteCombinations_Correct()
    {
        var calculator = new SampleSI16MagnitudeCalculator();
        for (int i = sbyte.MinValue; i <= sbyte.MaxValue; i++)
            for (int q = sbyte.MinValue; q <= sbyte.MaxValue; q++)
            {
                var sample = new SampleSI16((sbyte)i, (sbyte)q);
                var expected = (float)Math.Sqrt(i * i + q * q);

                var actual = calculator.GetMagnitude(sample);

                // Проверяем с небольшим допуском из-за float
                Assert.IsTrue(
                    condition: Math.Abs(expected - actual) < 1e-5f,
                    message: $"I={i}, Q={q}, expected={expected}, actual={actual}");
            }
    }
}
