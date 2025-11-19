using MathCore.DSP.Samples;

namespace MathCore.DSP.Tests.Samples;

/// <summary>Проверяет корректность вычисления аргумента для всех комбинаций I и Q типа sbyte</summary>
[TestClass]
public class SampleSI16ArgumentCalculatorTests
{
    /// <summary>Проверяет корректность вычисления аргумента для всех комбинаций I и Q типа sbyte</summary>
    [TestMethod]
    public void GetArgument_AllSByteCombinations_Correct()
    {
        var calculator = new SampleSI16ArgumentCalculator();
        for (int i = sbyte.MinValue; i <= sbyte.MaxValue; i++)
            for (int q = sbyte.MinValue; q <= sbyte.MaxValue; q++)
            {
                var sample = new SampleSI16((sbyte)i, (sbyte)q);
                var expected = (float)Math.Atan2(q, i);

                var actual = calculator.GetArgument(sample);

                // Проверяем с небольшим допуском из-за float
                Assert.IsLessThan(
1e-5f,
                    Math.Abs(expected - actual), $"I={i}, Q={q}, expected={expected}, actual={actual}");
            }
    }
}
