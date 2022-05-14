using System.ComponentModel;
using System.Diagnostics;

using MathCore.DSP.Infrastructure;

namespace MathCore.DSP.Filters;

/// <summary>Фильтр Чебышева нижних частот</summary>
public class ChebyshevHighPass : ChebyshevFilter
{
    private static (double[] A, double[] B) InitializeI(Specification Spec)
    {
        var N = (int)Math.Ceiling(arcch(Spec.kEps) / arcch(Spec.kW)); // Порядок фильтра
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var poles = GetNormedPolesI(N, Spec.EpsP);

        var high_pass_poles = TransformToHighPassW(poles, Spec.Ws);

        var z_poles = ToZArray(high_pass_poles, Spec.dt);

        var k_poles = z_poles.Multiply(z => (1 + z) / 2).Re;

        var g_norm = N.IsEven()
            ? Spec.Gp * k_poles
            : 1 * k_poles;

        var B = new double[N + 1];
        for (var i = 0; i < B.Length; i++)
            B[i] = SpecialFunctions.BinomialCoefficient(N, i) * (i % 2 == 0 ? g_norm : -g_norm);
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    private static (double[] A, double[] B) InitializeII(Specification Spec)
    {
        var N = (int)Math.Ceiling(arcch(Spec.kEps) / arcch(Spec.kW)); // Порядок фильтра
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");
        var (zeros, poles) = GetNormedPolesII(N, Spec.EpsS);

        var high_pass_zeros = TransformToHighPassW(zeros, Spec.Ws);
        if (N.IsOdd())
            high_pass_zeros = high_pass_zeros.Prepend(0);
        var high_pass_poles = TransformToHighPassW(poles, Spec.Ws);

        var z_zeros = ToZArray(high_pass_zeros, Spec.dt);
        var z_poles = ToZArray(high_pass_poles, Spec.dt);

        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToRe();
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        var k_zeros = z_zeros.Multiply(z => 1 + z).Re;
        var k_poles = z_poles.Multiply(z => 1 + z).Re;

        var g_norm = N.IsEven()
            ? Spec.Gp * k_poles / k_zeros
            : 1 * k_poles / k_zeros;

        B.Multiply(g_norm);

        return (A, B);
    }

    private static (double[] A, double[] B) InitializeIICorrected(Specification Spec)
    {
        var N = (int)Math.Ceiling(arcch(Spec.kEps) / arcch(Spec.kW)); // Порядок фильтра
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var (zeros, poles) = GetNormedPolesII(N, Spec.EpsS, Spec.kw);

        var high_pass_zeros = TransformToHighPassW(zeros, Spec.Ws);
        if (N.IsOdd())
            high_pass_zeros = high_pass_zeros.Prepend(0);
        var high_pass_poles = TransformToHighPassW(poles, Spec.Ws);

        var z_zeros = ToZArray(high_pass_zeros, Spec.dt);
        var z_poles = ToZArray(high_pass_poles, Spec.dt);

        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToRe();
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        var k_zeros = z_zeros.Multiply(z => 1 + z).Re;
        var k_poles = z_poles.Multiply(z => 1 + z).Re;

        var g_norm = N.IsEven()
            ? Spec.Gp * k_poles / k_zeros
            : 1 * k_poles / k_zeros;

        B.Multiply(g_norm);

        return (A, B);
    }

    /// <summary>Инициализация нового фильтра Чебышева нижних частот</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="fs">Частота заграждения</param>
    /// <param name="fp">Частота пропускания</param>
    /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
    /// <param name="Gs">Затухание в полосе заграждения (0.01        = -40 дБ)</param>
    /// <param name="Type">Тип (род) фильтра чебышева</param>
    public ChebyshevHighPass(double dt,
        double fs,
        double fp,
        double Gp = 0.891250938,
        double Gs = 0.01,
        ChebyshevType Type = ChebyshevType.I)
        : this(GetSpecification(dt, fs, fp, Gp, Gs), Type) { }

    public ChebyshevHighPass(Specification Spec, ChebyshevType Type = ChebyshevType.I)
        : this(Type switch
        {
            ChebyshevType.I => InitializeI(Spec),
            ChebyshevType.II => InitializeII(Spec),
            ChebyshevType.IICorrected => InitializeIICorrected(Spec),
            _ => throw new InvalidEnumArgumentException(nameof(Type), (int)Type, typeof(ChebyshevType))
        }, Spec, Type) { }

    private ChebyshevHighPass((double[] A, double[] B) config, Specification Spec, ChebyshevType Type) : base(config.B, config.A, Spec, Type) { }
}