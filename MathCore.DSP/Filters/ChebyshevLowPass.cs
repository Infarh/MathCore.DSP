using System.ComponentModel;

using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions;

namespace MathCore.DSP.Filters;

/// <summary>Фильтр Чебышева нижних частот</summary>
public class ChebyshevLowPass : ChebyshevFilter
{
    private static (double[] A, double[] B) InitializeI(Specification Spec)
    {
        var N = (int)Ceiling(arcch(Spec.kEps) / arcch(Spec.kW)); // Порядок фильтра

        var poles = GetNormedPolesI(N, Spec.EpsP);
        var z_poles = ToZArray(poles, Spec.dt, Spec.Wp);

        var A = GetCoefficientsInverted(z_poles).ToRe();

        var g_norm = (N.IsOdd() ? 1 : Spec.Gp)
            / (2.Power(N) / z_poles.Multiply(z => 1 - z).Re);

        var B = Enumerable
           .Range(0, N + 1)
           .ToArray(i => g_norm * BinomialCoefficient(N, i));

        return (A!, B!);
    }

    private static (double[] A, double[] B) InitializeII(Specification Spec)
    {
        var N = (int)Ceiling(arcch(Spec.kEps) / arcch(Spec.kW)); // Порядок фильтра
        var (zeros, poles) = GetNormedPolesII(N, Spec.EpsS);

        var z_zeros = N.IsEven()
            ? ToZArray(zeros, Spec.dt, Spec.Wp)
            : ToZ(zeros, Spec.dt, Spec.Wp).Prepend(-1).ToArray();
        var z_poles = ToZArray(poles, Spec.dt, Spec.Wp);

        var B = GetCoefficientsInverted(z_zeros).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var g_norm = 1 / (z_zeros.Multiply(z => 1 - z).Re / z_poles.Multiply(z => 1 - z).Re);

        for (var i = 0; i < B!.Length; i++)
            B[i] *= g_norm;

        return (A, B);
    }

    private static (double[] A, double[] B) InitializeIICorrected(Specification Spec)
    {
        var N = (int)Ceiling(arcch(Spec.kEps) / arcch(Spec.kW)); // Порядок фильтра
        var (zeros, poles) = GetNormedPolesII(N, Spec.EpsS);

        var kw = Spec.kw;
        var z_zeros = N.IsEven()
            ? ToZArray(zeros, Spec.dt, Spec.Wp * kw)
            : ToZ(zeros, Spec.dt, Spec.Wp * kw).Prepend(-1).ToArray();
        var z_poles = ToZArray(poles, Spec.dt, Spec.Wp * kw);

        var B = GetCoefficientsInverted(z_zeros).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var g_norm = 1 / (z_zeros.Multiply(z => 1 - z).Re / z_poles.Multiply(z => 1 - z).Re);

        for (var i = 0; i < B.Length; i++)
            B[i] *= g_norm;

        return (A, B);
    }

    /// <summary>Инициализация нового фильтра Чебышева нижних частот</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="fp">Частота пропускания</param>
    /// <param name="fs">Частота заграждения</param>
    /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
    /// <param name="Gs">Затухание в полосе заграждения (0.005623413 = -45 дБ)</param>
    /// <param name="Type">Тип (род) фильтра чебышева</param>
    public ChebyshevLowPass(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.005623413, ChebyshevType Type = ChebyshevType.I)
        : this(GetSpecification(dt, fp, fs, Gp, Gs), Type) { }

    public ChebyshevLowPass(Specification Spec, ChebyshevType Type = ChebyshevType.I)
        : this(Type switch
        {
            ChebyshevType.I => InitializeI(Spec),
            ChebyshevType.II => InitializeII(Spec),
            ChebyshevType.IICorrected => InitializeIICorrected(Spec),
            _ => throw new InvalidEnumArgumentException(nameof(Type), (int)Type, typeof(ChebyshevType))
        }, Spec, Type) { }

    private ChebyshevLowPass((double[] A, double[] B) config, Specification Spec, ChebyshevType Type) : base(config.B, config.A, Spec, Type) { }
}