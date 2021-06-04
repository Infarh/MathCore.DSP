// ReSharper disable once CheckNamespace
namespace System
{
    internal static class IntExtensions
    {
        public static (int div, int mod) DivMod(this int value, int n) => (value / n, value % n);
    }
}
