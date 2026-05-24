using MathCore.DSP.Filters;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ComplexIIRTests
{
    [TestMethod]
    public void ComplexIIR_Process_MatchesIndependentChannelFiltering()
    {
        var random = new Random(103);
        var samples = Enumerable
           .Range(0, 256)
           .Select(_ => new Complex(random.NextDouble() * 2 - 1, random.NextDouble() * 2 - 1))
           .ToArray();

        var a = new[] { 1.0, -0.81, 0.17 };
        var b = new[] { 0.14, 0.05, 0.02 };

        var expected_i = samples.Select(s => s.Re).FilterIIR(a, b).ToArray();
        var expected_q = samples.Select(s => s.Im).FilterIIR(a, b).ToArray();

        var filter = new ComplexIIR(b, a);
        var actual = filter.Process(samples).ToArray();

        Assert.AreEqual(samples.Length, actual.Length);

        const double eps = 1e-12;
        for (var i = 0; i < actual.Length; i++)
        {
            Assert.AreEqual(expected_i[i], actual[i].Re, eps, $"Неверная I-компонента на позиции {i}");
            Assert.AreEqual(expected_q[i], actual[i].Im, eps, $"Неверная Q-компонента на позиции {i}");
        }
    }

    [TestMethod]
    public void ComplexIIR_Reset_RestoresInitialState()
    {
        var a = new[] { 1.0, -0.72, 0.13 };
        var b = new[] { 0.21, 0.08, 0.03 };
        var filter = new ComplexIIR(b, a);

        var input = new[]
        {
            new Complex(0.25, -0.5),
            new Complex(-1.2, 0.75),
            new Complex(0.4, 0.1),
            new Complex(0, 1)
        };

        var first = filter.Process(input).ToArray();
        filter.Reset();
        var second = filter.Process(input).ToArray();

        CollectionAssert.AreEqual(first, second);
    }
}
