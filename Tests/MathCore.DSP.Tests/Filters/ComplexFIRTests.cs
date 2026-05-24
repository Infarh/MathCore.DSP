using MathCore.DSP.Filters;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class ComplexFIRTests
{
    [TestMethod]
    public void ComplexFIR_Process_MatchesIndependentChannelFiltering()
    {
        var random = new Random(101);
        var samples = Enumerable
           .Range(0, 256)
           .Select(_ => new Complex(random.NextDouble() * 2 - 1, random.NextDouble() * 2 - 1))
           .ToArray();

        double[] impulse_response = [0.15, 0.35, 0.25, 0.1];

        var expected_i = samples.Select(s => s.Re).FilterFIR(impulse_response).ToArray();
        var expected_q = samples.Select(s => s.Im).FilterFIR(impulse_response).ToArray();

        var filter = new ComplexFIR(impulse_response);
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
    public void ComplexFIR_Reset_RestoresInitialState()
    {
        double[] impulse_response = [0.3, 0.2, 0.1];
        var filter = new ComplexFIR(impulse_response);

        var input = new[]
        {
            new Complex(1, -2),
            new Complex(0.5, 3),
            new Complex(-1, 0.25)
        };

        var first = filter.Process(input).ToArray();
        filter.Reset();
        var second = filter.Process(input).ToArray();

        CollectionAssert.AreEqual(first, second);
    }
}
