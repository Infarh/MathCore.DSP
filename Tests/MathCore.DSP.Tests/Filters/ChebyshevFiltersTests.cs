using Microsoft.VisualStudio.TestTools.UnitTesting;

using static System.Math;

namespace MathCore.DSP.Tests.Filters;

public abstract class ChebyshevFiltersTests : UnitTest
{
    protected static double arcsh(double x) => Log(x + Sqrt(x * x + 1));
    protected static double arcch(double x) => Log(x + Sqrt(x * x - 1));
}
