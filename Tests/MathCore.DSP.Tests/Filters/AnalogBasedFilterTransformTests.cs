using MathCore.DSP.Filters;
using MathCore.DSP.Tests.Infrastructure;

namespace MathCore.DSP.Tests.Filters;

[TestClass, TestCategory("123")]
public class AnalogBasedFilterTransformTests : UnitTest
{
    [TestMethod]
    public void LowToLow_transform_poles()
    {
        Complex[] poles =
        [
            -0.364,
            (-0.057, 0.997),
            (-0.224, 0.716)
        ];

        const double wp = 2 * Math.PI;

        var transformed_poles = AnalogBasedFilter.Transform.ToLow(poles, wp).ToArray();

        transformed_poles[0].AssertEquals(poles[0] * wp);
        transformed_poles[1].AssertEquals(poles[1] * wp);
        transformed_poles[2].AssertEquals(poles[2] * wp);
    }

    [TestMethod]
    public void LowToHigh_transform_poles()
    {
        Complex[] poles =
        [
            -0.364,
            (-0.057, 0.997),
            (-0.224, 0.716)
        ];

        const double wp = 2 * Math.PI;

        var transformed_poles = AnalogBasedFilter.Transform.ToHigh(poles, wp).ToArray();

        transformed_poles[0].AssertEquals(wp / poles[0]);
        transformed_poles[1].AssertEquals(wp / poles[1]);
        transformed_poles[2].AssertEquals(wp / poles[2]);
    }

    [TestMethod]
    public void LowToBandPass_transform_poles()
    {
        Complex[] poles =
        [
            -0.364,
            (-0.057, 0.997),
            (-0.224, 0.716)
        ];

        var zeros = Array.Empty<Complex>();

        const double fpl = 2 / Consts.pi2;
        const double fph = 12 / Consts.pi2;

        var (transformed_poles, _) = AnalogBasedFilter.Transform.ToBandPass(poles, zeros, fpl, fph);
        var expected_poles = AnalogBasedFilter.TransformToBandPassW(poles, Consts.pi2 * fpl, Consts.pi2 * fph).ToArray();

        transformed_poles.AssertEquals(AccuracyComplex.Eps(1e-12), expected_poles);
    }

    [TestMethod]
    public void LowToBandStop_transform_poles()
    {
        Complex[] poles =
        [
            -0.364,
            (-0.057, 0.997),
            (-0.224, 0.716)
        ];

        var zeros = Array.Empty<Complex>();

        const double fpl = 2 / Consts.pi2;
        const double fph = 12 / Consts.pi2;

        var (transformed_poles, _) = AnalogBasedFilter.Transform.ToBandStop(poles, zeros, fpl, fph);
        var expected_poles = AnalogBasedFilter.TransformToBandStopW(poles, Consts.pi2 * fpl, Consts.pi2 * fph).ToArray();

        transformed_poles.AssertEquals(AccuracyComplex.Eps(1e-12), expected_poles);
    }
}