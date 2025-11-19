#if !NET8_0_OR_GREATER
using System.Runtime.CompilerServices;

namespace MathCore.DSP.Infrastructure;

internal static class ArgumentNullExceptionEx
{
    extension(ArgumentNullException)
    {
        public static void ThrowIfNull(object? argument, [CallerArgumentExpression("argument")] string? ParamName = null)
        {
            if (argument != null)
                return;
            ArgumentNullException.Throw(ParamName);
        }

        [DoesNotReturn]
        public static void Throw(string ParamName) => throw new ArgumentNullException(ParamName);
    }
}

[AttributeUsage(AttributeTargets.Method, Inherited = false)]
public sealed class DoesNotReturnAttribute : Attribute;

#endif