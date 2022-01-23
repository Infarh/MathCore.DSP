using System;
using System.Linq;

using MathCore.DSP.Filters;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Filters;

[TestClass, TestCategory("123")]
public class AnalogBasedFilterTransformTests : UnitTest
{
    [TestMethod]
    public void LowToLow_transform_poles()
    {
        Complex[] poles =
        {
            -0.364,
            (-0.057, 0.997),
            (-0.224, 0.716)
        };

        var poles_lowpas = AnalogBasedFilter.Transform.ToLow(poles, 2 * Math.PI).ToArray();
    }
}