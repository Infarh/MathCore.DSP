namespace MathCore.DSP.Samples;

public abstract class SampleSI16Calculator
{
    protected static int GetIndex(int I, int Q)
    {
        if (I < Q) // Используем симметрию: меняем I и Q
        {
            // ReSharper disable once SwapViaDeconstruction
            var tmp = I;
            I = Q;
            Q = tmp;
        }

        return (I * (I + 1)) / 2 + Q; // Вычисляем индекс для I >= Q
    }
}