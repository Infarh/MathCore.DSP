// ReSharper disable CheckNamespace
namespace System.Runtime.CompilerServices;

[AttributeUsage(AttributeTargets.Parameter)]
internal sealed class CallerArgumentExpressionAttribute(string ParameterName) : Attribute
{
    public string ParameterName { get; } = ParameterName;
}