using System.Diagnostics;

using MathCore.DSP.Infrastructure;

using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions;

namespace MathCore.DSP.Filters;

public class ButterworthHighPass : ButterworthFilter
{
    public static (double[] A, double[] B) GetPolynoms(Specification Spec)
    {
        // Порядок фильтра
        var N = (int)Ceiling(Log(Spec.kEps) / -Log(Spec.kW));
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");
        var poles = GetNormPoles(N, Spec.EpsP).ToArray();
        poles.ToDebugEnum();

        // Масштабируем полюса на требуемую частоту пропускания
        var high_pass_poles = TransformToHighPassW(poles, Spec.Wp);

        // Переходим из p-плоскости в z-плоскость
        var z_poles = ToZArray(high_pass_poles, Spec.dt);

        // Вычисляем нормирующий множитель
        var g_norm = z_poles.Multiply(z => (1 + z) / 2).Abs;
        g_norm.ToDebug();

        var B = new double[N + 1];
        for (var i = 0; i < B.Length; i++)
            B[i] = BinomialCoefficient(N, i) * (i % 2 == 0 ? g_norm : -g_norm);
        var A = GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    public ButterworthHighPass(
        double dt,
        double fs,
        double fp,
        double Gp = 0.891250938,
        double Gs = 0.031622777)
        : this(GetSpecification(dt, fp, fs, Gp, Gs)) { }

    public ButterworthHighPass(Specification Spec) : this(GetPolynoms(Spec), Spec) { }

    public ButterworthHighPass((double[] A, double[] B) config, Specification Spec) : base(config.B, config.A, Spec) { }
}