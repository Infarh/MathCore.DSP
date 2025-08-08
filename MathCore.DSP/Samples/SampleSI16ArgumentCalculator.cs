namespace MathCore.DSP.Samples;

public class SampleSI16ArgumentCalculator : SampleSI16Calculator
{
    private readonly float[] _Argument = new float[8385];

    public SampleSI16ArgumentCalculator()
    {
        var index = 0;
        for (var i = 0; i <= 128; i++)
            for (var q = 0; q <= i; q++)
            {
#if NET8_0_OR_GREATER
                float re = i;
                float im = q;
                _Argument[index] = i == 0 && q == 0
                    ? 0
                    : MathF.Atan2(im, re);
#else
            double re = i;
            double im = q;
            _Argument[index] = i == 0 && q == 0
                ? 0
                : (float)Math.Atan2(im, re);
#endif
                index++;
            }
    }

    public float GetArgument(SampleSI16 Sample)
    {
        var (re, im) = Sample;

        if (re == 0 && im == 0)
            return 0;

        var abs_i = Math.Abs((int)re);
        var abs_q = Math.Abs((int)im);

        var swapped = abs_i < abs_q;
        var index = GetIndex(abs_i, abs_q);
        var arg = _Argument[index];

#if NET8_0_OR_GREATER
        if (swapped) // Если I < Q, то аргумент для (Q, I) = π/2 - arg для (I, Q)
            arg = MathF.PI / 2 - arg;

        arg = (re, im) switch
        {
            ( >= 0, >= 0) => arg,           // I квадрант: без изменений
            ( < 0, >= 0) => MathF.PI - arg, // II квадрант: π - arg
            ( < 0, < 0) => arg - MathF.PI,  // III квадрант: arg - π
            _ => -arg                       // IV квадрант: -arg
        };

        // Нормализация аргумента в диапазон (-π, π]
        var correction = arg switch
        {
            <= -MathF.PI => 2 * MathF.PI,   // Если arg <= -π, добавляем 2π
            > MathF.PI => -(2 * MathF.PI), // Если arg > π, вычитаем 2π
            _ => 0
        };
        return arg + correction;
#else
        if (swapped) // Если I < Q, то аргумент для (Q, I) = π/2 - arg для (I, Q)
            arg = (float)(Math.PI / 2 - arg);

        arg = (re, im) switch
        {
            ( >= 0, >= 0) => arg,                   // I квадрант: без изменений
            ( < 0, >= 0) => (float)(Math.PI - arg), // II квадрант: π - arg
            ( < 0, < 0) => (float)(arg - Math.PI),  // III квадрант: arg - π
            _ => -arg                               // IV квадрант: -arg
        };

        // Нормализация аргумента в диапазон (-π, π]
        return arg + arg switch
        {
            <= -(float)Math.PI => (float)(2 * Math.PI), // Если arg <= -π, добавляем 2π
            > (float)Math.PI => -(float)(2 * Math.PI), // Если arg > π, вычитаем 2π
            _ => 0
        };
#endif
    }
}