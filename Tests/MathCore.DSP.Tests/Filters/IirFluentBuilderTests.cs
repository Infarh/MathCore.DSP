using MathCore.DSP.Filters;
using MathCore.DSP.Filters.Builders;

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

    [TestMethod]
    public void CreateBySamplingFrequencyFactory()
    {
        const double fd = 5_000;

        var filter = Filter.IIRBySamplingFrequency(fd)
            .Butterworth
            .LowPass(500, 1_500)
            .Create();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.ButterworthLowPass));
    }

    [TestMethod]
    public void CreateBandPassByCenterWithDbGains()
    {
        const double dt = 1d / 10_000;

        var filter = Filter.IIR(dt)
            .Elliptic
            .BandPassByCenter(1_500, 1_000, 500)
            .WithGainsInDb(1, 40)
            .Create();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.EllipticBandPass));
    }

    [TestMethod]
    public void TryCreateReturnsFalseForInvalidSpecification()
    {
        const double dt = 1d / 10_000;

        var is_created = Filter.IIR(dt)
            .Butterworth
            .LowPass(2_000, 1_000)
            .TryCreate(out var filter);

        Assert.IsFalse(is_created);
        Assert.IsNull(filter);
    }

    [TestMethod]
    public void CreateBandStopByCenter()
    {
        const double dt = 1d / 10_000;

        var filter = Filter.IIR(dt)
            .Chebyshev(ChebyshevType.II)
            .BandStopByCenter(1_500, 1_000, 500)
            .Create();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.ChebyshevBandStop));
    }

    [TestMethod]
    public void CreateByDslWithSamplingFrequencyPassStopBands()
    {
        var filter = Filter.IIR()
            .WithSamplingFrequency(10_000)
            .Butterworth
            .BandPass()
            .WithPassband(1_000, 2_000)
            .WithStopband(500, 2_500)
            .ToSpecification()
            .Create();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.ButterworthBandPass));
    }

    [TestMethod]
    public void CreateByDslWithSamplingAndDbConfiguration()
    {
        var filter = Filter.IIR()
            .WithSampling(1d / 10_000)
            .Elliptic
            .LowPass()
            .WithPassband(1_000)
            .WithStopband(2_000)
            .ToSpecification()
            .WithGainsInDb(1, 40)
            .Create();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.EllipticLowPass));
    }

    [TestMethod]
    public void GetSpecificationReturnsButterworthLowPassSpecification()
    {
        const double dt = 1d / 5_000;

        var specification = Filter.IIR(dt)
            .Butterworth
            .LowPass(500, 1_500)
            .GetSpecification();

        Assert.IsInstanceOfType(specification, typeof(ButterworthLowPassSpecification));
    }

    [TestMethod]
    public void SpecificationBuildsFilter()
    {
        const double dt = 1d / 10_000;

        var specification = Filter.IIR(dt)
            .Chebyshev(ChebyshevType.II)
            .BandStop(500, 1_000, 2_000, 2_500)
            .GetSpecification();

        var filter = specification.CreateFilter();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.ChebyshevBandStop));
    }

    [TestMethod]
    public void DslSpecificationImplementsHorizontalInterfaces()
    {
        var specification = Filter.IIR()
            .WithSamplingFrequency(10_000)
            .Elliptic
            .BandPass()
            .WithPassband(1_000, 2_000)
            .WithStopband(500, 2_500)
            .GetSpecification();

        Assert.IsInstanceOfType(specification, typeof(IIirBandFrequenciesSpecification));
        Assert.IsInstanceOfType(specification, typeof(IIirGainSpecification));
        Assert.IsInstanceOfType(specification, typeof(IIirFamilySpecification));
    }

    [TestMethod]
    public void RcBuilderReturnsSpecificationAndBuildsFilter()
    {
        const double dt = 1d / 1_000;

        var specification = Filter.IIR(dt)
            .RC
            .HighPass(50)
            .GetSpecification();

        var filter = specification.CreateFilter();

        Assert.IsInstanceOfType(specification, typeof(RcHighPassSpecification));
        Assert.IsInstanceOfType(filter, typeof(RCHighPass));
    }

    [TestMethod]
    public void BuildSpecificationAliasReturnsSpecification()
    {
        const double dt = 1d / 5_000;

        var specification = Filter.IIR(dt)
            .Butterworth
            .LowPass(500, 1_500)
            .BuildSpecification();

        Assert.IsInstanceOfType(specification, typeof(ButterworthLowPassSpecification));
    }

    [TestMethod]
    public void BuildFilterAliasReturnsFilter()
    {
        var filter = Filter.IIR()
            .WithSamplingFrequency(10_000)
            .Elliptic
            .BandPass()
            .WithPassband(1_000, 2_000)
            .WithStopband(500, 2_500)
            .BuildFilter();

        Assert.IsInstanceOfType(filter, typeof(global::MathCore.DSP.Filters.EllipticBandPass));
    }

    [TestMethod]
    public void TryGetSpecificationReturnsFalseForIncompleteDsl()
    {
        var success = Filter.IIR()
            .WithSamplingFrequency(10_000)
            .Butterworth
            .LowPass()
            .WithPassband(1_000)
            .TryGetSpecification(out var specification);

        Assert.IsFalse(success);
        Assert.IsNull(specification);
    }
}
