namespace MathCore.DSP.Signals.Interpolators;

/// <summary>Линейный интерполятор</summary>
[Serializable]
public class Linear : Interpolator
{
    public override double GetValue(double[] Values, double pos)
    {
        var N = Values.Length;
        var x1 = (int)pos;
        if(x1 >= N - 1) 
            return Values[N - 1];

        var x2 = x1 + 1;
        if(x2 <= 0) 
            return Values[0];

        var y1 = Values[x1];
        var y2 = Values[x2];

        return (y2 - y1) * (pos - x1) / (x2 - x1) + y1;
    }

    public override double SetValue(double[] Values, double pos, double Value) => throw new NotSupportedException();
}