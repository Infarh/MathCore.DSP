// ReSharper disable InconsistentNaming

using System.Runtime.InteropServices;
// ReSharper disable ArrangeNullCheckingPattern

namespace MathCore.DSP.Samples;

/// <summary>Значение отсчёта квадратурного цифрового сигнала с целочисленными 8-битовыми компонентами</summary>
/// <param name="I">Синфазная компонента</param>
/// <param name="Q">Квадратурная компонента</param>
[StructLayout(LayoutKind.Sequential)]
public readonly record struct SampleSI16(sbyte I, sbyte Q)
{
#if NET8_0_OR_GREATER
    /// <summary>Функция вычисления модуля</summary>
    public static Func<SampleSI16, float> GetAbs
    {
        get;
        set
        {
            if (value is not { })
                value = s => MathF.Sqrt(s.I * s.I + s.Q * s.Q);
            field = value;
        }
    } = s => MathF.Sqrt(s.I * s.I + s.Q * s.Q);

    /// <summary>Функция вычисления аргумента</summary>
    public static Func<SampleSI16, float> GetArg
    {
        get;
        set
        {
            if (value is not { })
                value = s => MathF.Atan2(s.Q, s.I);
            field = value;
        }
    } = s => MathF.Atan2(s.Q, s.I);
#else
    /// <summary>Функция вычисления модуля</summary>
    public static Func<SampleSI16, float> GetAbs
    {
        get;
        set
        {
            if (value is not { })
                value = s => (float)Math.Sqrt(s.I * s.I + s.Q * s.Q);
            field = value;
        }
    } = s => (float)Math.Sqrt(s.I * s.I + s.Q * s.Q);

    /// <summary>Функция вычисления аргумента</summary>
    public static Func<SampleSI16, float> GetArg
    {
        get;
        set
        {
            if (value is not { })
                value = s => (float)Math.Atan2(s.Q, s.I);
            field = value;
        }
    } = s => (float)Math.Atan2(s.Q, s.I);
#endif

    /// <summary>Модуль</summary>
    public float Abs => GetAbs(this);

    /// <summary>Аргумент (фаза)</summary>
    public float Arg => GetArg(this);
}
