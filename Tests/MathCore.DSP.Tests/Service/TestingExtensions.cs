using System.Collections.Generic;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Service
{
    internal static class TestingExtensions
    {
        public static ValueChecker<T> AssertThanValue<T>(this T value) => Assert.That.Value(value);

        public static CollectionChecker<T> AssertThatCollection<T>(this ICollection<T> collection) => Assert.That.Collection(collection);
    }
}
