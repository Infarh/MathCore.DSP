using MathCore.DSP.Filters;

namespace MathCore.DSP.Tests.Filters;

[TestClass]
public class IirFluentBuilderTests : UnitTest
{
    [TestMethod]
    public void CreateButterworthLowPassByFluentBuilder()
    {
        const double dt = 1d / 5_000;

        var filter = Filter.IIR(dt)
            .Butterworth
            .LowPass(500, 1_500)
            .Create();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.ButterworthLowPass));
    }

    [TestMethod]
    public void CreateChebyshevBandPassByFluentBuilderWithConfiguredType()
    {
        const double dt = 1d / 10_000;

        var filter = Filter.IIR(dt)
            .Chebyshev()
            .BandPass(500, 1_000, 2_000, 2_500)
            .WithChebyshevType(ChebyshevType.II)
            .Create();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.ChebyshevBandPass));
    }

    [TestMethod]
    public void CreateEllipticBandStopByFluentBuilder()
    {
        const double dt = 1d / 10_000;

        var filter = Filter.IIR(dt)
            .Elliptic
            .BandStop(500, 1_000, 2_000, 2_500)
            .WithGains(0.891250938, 0.01)
            .Create();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.EllipticBandStop));
    }

    [TestMethod]
    public void CreateRcAndRlcFiltersByFluentBuilder()
    {
        const double dt = 1d / 1_000;

        var rc_high_pass = Filter.IIR(dt).RC.HighPass(50).Create();
        var rlc_band_pass = Filter.IIR(dt).RLC.BandPass(120, 20).Create();

        Assert.IsInstanceOfType<RCHighPass>(rc_high_pass);
        Assert.IsInstanceOfType<RLCBandPass>(rlc_band_pass);
    }

    [TestMethod]
    public void CreateManyReturnsIndependentInstances()
    {
        const double dt = 1d / 5_000;

        var filters = Filter.IIR(dt)
            .Butterworth
            .HighPass(500, 1_000)
            .CreateMany(2)
            .ToArray();

        Assert.AreEqual(2, filters.Length);
        Assert.AreNotSame(filters[0], filters[1]);
    }
}
