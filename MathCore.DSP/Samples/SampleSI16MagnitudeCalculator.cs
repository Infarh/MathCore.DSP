namespace MathCore.DSP.Samples;

public class SampleSI16MagnitudeCalculator : SampleSI16Calculator
{
    private readonly float[] _Magnitude = new float[8385];

    public SampleSI16MagnitudeCalculator()
    {
        var index = 0;
        for (var i = 0; i <= 128; i++)
            for (var q = 0; q <= i; q++)
            {
#if NET8_0_OR_GREATER
                float re = i;
                float im = q;

                _Magnitude[index] = MathF.Sqrt(re * re + im * im);
#else
                double re = i;
                double im = q;

                _Magnitude[index] = (float)Math.Sqrt(re * re + im * im);
#endif
                index++;
            }
    }

    public float GetMagnitude(SampleSI16 s)
    {
        var (i, q) = s;
        var abs_i = Math.Abs((int)i);
        var abs_q = Math.Abs((int)q);
        var index = GetIndex(abs_i, abs_q);
        return _Magnitude[index];
    }
}