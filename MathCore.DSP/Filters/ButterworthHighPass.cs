using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions;

namespace MathCore.DSP.Filters;

public class ButterworthHighPass : ButterworthFilter
{
    public static (double[] A, double[] B) GetPolynoms(Specification Spec)
    {
        // Порядок фильтра
        var N = (int)Ceiling(Log(Spec.kEps) / Log(Spec.kW));
        var poles = GetNormPoles(N, Spec.EpsP).ToArray();

        // Масштабируем полюса на требуемую частоту пропускания
        var Wp = Spec.Wp;
        var translated_poles = TransformToHighPass(poles, Wp);
        var zeros = Enumerable.Repeat(Complex.ReValue(-1), N);

        // Переходим из p-плоскости в z-плоскость
        var dt = Spec.dt;
        var z_poles = ToZArray(translated_poles, dt);   
        var z_zeros = ToZArray(zeros, dt);
        //var z_poles = translated_poles.ToArray(p => ToZ(p, dt));
        // Вычисляем нормирующий множитель
        //var kz = GetNormalizeCoefficient(translated_poles, dt);
        var kz = poles.Multiply();

        var WpN = Wp.Pow(N);
        var k = WpN * kz / Spec.EpsP;

        var B = GetCoefficientsInverted(z_zeros).ToArray(v => (v * k).Re);
        var A = GetCoefficientsInverted(z_poles).ToRe()!;

        return (A, B);
    }

    public ButterworthHighPass(
        double dt,
        double fp,
        double fs,
        double Gp = 0.891250938,
        double Gs = 0.031622777) 
        : this(GetSpecification(dt, fp, fs, Gp, Gs)) { }

    public ButterworthHighPass(Specification Spec) : this(GetPolynoms(Spec), Spec) { }

    public ButterworthHighPass((double[] A, double[] B) config, Specification Spec) : base(config.B, config.A, Spec) { }
}