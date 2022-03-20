using System.Collections.Generic;

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathCore.DSP.Tests.Service;

internal static class TestingExtensions
{
    public static ValueChecker<T> AssertThatValue<T>(this T value) => Assert.That.Value(value);
    public static ValueChecker<T> AssertEquals<T>(this T value, T ActualValue) => Assert.That.Value(value).IsEqual(ActualValue);
    public static ValueChecker<double> AssertEquals(this double value, double ActualValue, double Eps) => Assert.That.Value(value).IsEqual(ActualValue, Eps);

    public static CollectionChecker<T> AssertThatCollection<T>(this ICollection<T> collection) => Assert.That.Collection(collection);

    public static CollectionChecker<T> AssertEquals<T>(this ICollection<T> collection, params T[] args) =>
        Assert.That.Collection(collection).IsEqualTo(args);

    public static EnumerableChecker<T> AssertThatEnumerable<T>(this IEnumerable<T> items) => Assert.That.Enumerable(items);

    public static EnumerableChecker<T> AssertEquals<T>(this IEnumerable<T> items, params T[] values) => Assert.That.Enumerable(items).IsEqualTo(values);
}