using System.ComponentModel;
using System.Diagnostics;

using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters;

/// <summary>Фильтр Чебышева нижних частот</summary>
public class ChebyshevLowPass : ChebyshevFilter
{
    public static double GetFrequencyStopTypeI(double dt, double fp, int Order, double Gp = 0.891250938, double Gs = 0.01)
    {
        var g2 = (1 / (Gs * Gs) - 1) / (1 / (Gp * Gp) - 1);
        //var pi_dt = PI * dt;
        var Fp = ToAnalogFrequency(fp, dt);
        //var kW = Tan(fs * pi_dt) / Tan(fp * pi_dt);
        var log_g = Math.Log(Sqrt(g2) + Sqrt(g2 - 1));
        return Fp * Cosh(log_g / Order);
    }
    public static double GetFrequencyPassTypeI(double dt, double fs, int Order, double Gp = 0.891250938, double Gs = 0.01)
    {
        var g2 = (1 / (Gs * Gs) - 1) / (1 / (Gp * Gp) - 1);
        //var pi_dt = PI * dt;
        var Fs = ToAnalogFrequency(fs, dt);
        //var kW = Tan(fs * pi_dt) / Tan(fp * pi_dt);
        var log_g = Math.Log(Sqrt(g2) + Sqrt(g2 - 1));
        return Fs / Cosh(log_g / Order);
    }

    public static double GetGsTypeI(double dt, double fp, double fs, int Order, double Gp = 0.891250938)
    {
        var pi_dt = PI * dt;
        var f = Tan(fs * pi_dt) / Tan(fp * pi_dt);

        var ff = (f + Sqrt(f * f - 1)).Pow(Order);

        var q = ff + 1 / ff;
        return 1 / Sqrt(0.25 * q * q * (1 / (Gp * Gp) - 1) + 1);
    }

    public static double GetGpTypeI(double dt, double fp, double fs, int Order, double Gs = 0.01)
    {
        var pi_dt = PI * dt;
        var f = Tan(fs * pi_dt) / Tan(fp * pi_dt);

        var ff = (f + Sqrt(f * f - 1)).Pow(Order);

        var q = ff + 1 / ff;
        return 1 / Sqrt((4 / (Gs * Gs) - 4) / (q * q));
    }

    public static int GetOrderTypeI(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.01)
    {
        var kEps = Sqrt((1 / (Gs * Gs) - 1) / (1 / (Gp * Gp) - 1));
        var pi_dt = PI * dt;
        var kW = Tan(fs * pi_dt) / Tan(fp * pi_dt);
        var N = (int)Ceiling(arcch(kEps) / arcch(kW));
        return N;
    }

    private static (double[] A, double[] B) InitializeI(Specification Spec)
    {
        var N = (int)Ceiling(arcch(Spec.kEps) / arcch(Spec.kW)); // Порядок фильтра
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var poles = GetNormedPolesI(N, Spec.EpsP);
        var z_poles = ToZArray(poles, Spec.dt, Spec.Wp);

        var A = GetCoefficientsInverted(z_poles).ToRe();

        var g_norm = (N.IsOdd() ? 1 : Spec.Gp)
            / (2.Power(N) / z_poles.Multiply(z => 1 - z).Re);

        var B = Enumerable
           .Range(0, N + 1)
           .ToArray(i => g_norm * BinomialCoefficient(N, i));

        return (A, B);
    }

    private static (double[] A, double[] B) InitializeII(Specification Spec)
    {
        var N = (int)Ceiling(arcch(Spec.kEps) / arcch(Spec.kW)); // Порядок фильтра
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");
        var (zeros, poles) = GetNormedPolesII(N, Spec.EpsS);

        var z_zeros = N.IsEven()
            ? ToZArray(zeros, Spec.dt, Spec.Wp)
            : ToZ(zeros, Spec.dt, Spec.Wp).Prepend(-1).ToArray();
        var z_poles = ToZArray(poles, Spec.dt, Spec.Wp);

        var B = GetCoefficientsInverted(z_zeros).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var g_norm = 1 / (z_zeros.Multiply(z => 1 - z).Re / z_poles.Multiply(z => 1 - z).Re);

        B.Multiply(g_norm);

        return (A, B);
    }

    private static (double[] A, double[] B) InitializeIICorrected(Specification Spec)
    {
        var N = (int)Ceiling(arcch(Spec.kEps) / arcch(Spec.kW)); // Порядок фильтра
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");
        var (zeros, poles) = GetNormedPolesII(N, Spec.EpsS);

        var kw = Spec.kw;
        var z_zeros = N.IsEven()
            ? ToZArray(zeros, Spec.dt, Spec.Wp * kw)
            : ToZ(zeros, Spec.dt, Spec.Wp * kw).Prepend(-1).ToArray();
        var z_poles = ToZArray(poles, Spec.dt, Spec.Wp * kw);

        var B = GetCoefficientsInverted(z_zeros).ToRe();
        var A = GetCoefficientsInverted(z_poles).ToRe();

        var g_norm = 1 / (z_zeros.Multiply(z => 1 - z).Re / z_poles.Multiply(z => 1 - z).Re);

        B.Multiply(g_norm);

        return (A, B);
    }

    /// <summary>Инициализация нового фильтра Чебышева нижних частот</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="fp">Частота пропускания</param>
    /// <param name="fs">Частота заграждения</param>
    /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
    /// <param name="Gs">Затухание в полосе заграждения (0.01        = -40 дБ)</param>
    /// <param name="Type">Тип (род) фильтра чебышева</param>
    public ChebyshevLowPass(
        double dt,
        double fp,
        double fs,
        double Gp = 0.891250938,
        double Gs = 0.01,
        ChebyshevType Type = ChebyshevType.I)
        : this(GetSpecification(dt, fp, fs, Gp, Gs), Type) { }

    public ChebyshevLowPass(Specification Spec, ChebyshevType Type = ChebyshevType.I)
        : this(Type switch
        {
            ChebyshevType.I => InitializeI(Spec),
            ChebyshevType.II => InitializeII(Spec),
            ChebyshevType.IICorrected => InitializeIICorrected(Spec),
            _ => throw new InvalidEnumArgumentException(nameof(Type), (int)Type, typeof(ChebyshevType))
        }, Spec, Type)
    { }

    private ChebyshevLowPass((double[] A, double[] B) config, Specification Spec, ChebyshevType Type) : base(config.B, config.A, Spec, Type) { }
}