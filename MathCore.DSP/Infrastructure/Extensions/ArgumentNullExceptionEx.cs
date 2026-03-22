#if NETSTANDARD2_0
using System.Runtime.CompilerServices;

namespace MathCore.DSP.Infrastructure.Extensions;

internal static class ArgumentNullExceptionEx
{
    extension(ArgumentNullException)
    {
        public static void ThrowIfNull(object? argument, [CallerArgumentExpression(nameof(argument))] string? ParamName = null)
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