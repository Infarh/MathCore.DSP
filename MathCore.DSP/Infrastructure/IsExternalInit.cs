#if !NET7_0_OR_GREATER

// ReSharper disable once CheckNamespace
namespace System.Runtime.CompilerServices;

/// <summary>
/// Dummy class so C# init-only properties can compile on NetStandard.
/// </summary>
internal sealed class IsExternalInit { }

#endif