using MathCore.DSP.Samples;

namespace MathCore.DSP.Tests.Samples;

[TestClass]
public class SampleSI16Tests
{
    [ClassInitialize]
    public static void Initialize(TestContext context)
    {
        var abs_calculator = new SampleSI16MagnitudeCalculator();
        var arg_calculator = new SampleSI16ArgumentCalculator();

        SampleSI16.GetAbs = abs_calculator.GetMagnitude;
        SampleSI16.GetArg = arg_calculator.GetArgument;
    }

    [TestMethod]
    public void AbsArgTest()
    {
        for (var i = -128; i <= 127; i++)
            for (var q = -128; q <= 127; q++)
            {
                var expected_abs = MathF.Sqrt(i * i + q * q);
                var expected_arg = i == 0 && q == 0 ? 0 : MathF.Atan2(q, i);

                var sample = new SampleSI16((sbyte)i, (sbyte)q);

                var actual_abs = sample.Abs;
                var actual_arg = sample.Arg;

                Assert.AreEqual(expected_abs, actual_abs);
                Assert.AreEqual(expected_arg, actual_arg, 2.5e-7);
            }
    }
}
