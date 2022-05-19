using System.Diagnostics;

using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions.EllipticJacobi;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters;

public class EllipticBandPass : EllipticFilter
{
    private static void CheckFrequencies(double dt, double fsl, double fpl, double fph, double fsh)
    {
        if (dt <= 0)
            throw new InvalidOperationException($"Период дискретизации dt={dt} не может быть меньше, либо равен нулю")
            {
                Data =
                {
                    { nameof(dt), dt },
                    { "fd", 1/dt },
                    { nameof(fsl), fsl },
                    { nameof(fpl), fpl },
                    { nameof(fph), fph },
                    { nameof(fsh), fsh },
                }
            };

        if (1 / dt == 0)
            throw new InvalidOperationException("Частота дискретизации не может быть равна нулю")
            {
                Data =
                {
                    { nameof(dt), dt },
                    { "fd", 1/dt },
                    { nameof(fsl), fsl },
                    { nameof(fpl), fpl },
                    { nameof(fph), fph },
                    { nameof(fsh), fsh },
                }
            };

        if (fsl >= fpl)
            throw new InvalidOperationException($"Нижняя частота среза fsl должна быть ниже нижней частоты пропускания fpl\r\n  dt={dt}\r\n  fsl={fsl}\r\n  fpl={fpl}\r\n  fph={fph}\r\n  fsh={fsh}")
            {
                Data =
                {
                    { nameof(dt), dt },
                    { "fd", 1/dt },
                    { nameof(fsl), fsl },
                    { nameof(fpl), fpl },
                    { nameof(fph), fph },
                    { nameof(fsh), fsh },
                }
            };

        if (fpl >= fph)
            throw new InvalidOperationException($"Нижняя частота пропускания fpl должна быть ниже верхней частоты пропускания fph\r\n  dt={dt}\r\n  fsl={fsl}\r\n  fpl={fpl}\r\n  fph={fph}\r\n  fsh={fsh}")
            {
                Data =
                {
                    { nameof(dt), dt },
                    { "fd", 1/dt },
                    { nameof(fsl), fsl },
                    { nameof(fpl), fpl },
                    { nameof(fph), fph },
                    { nameof(fsh), fsh },
                }
            };

        if (fph >= fsh)
            throw new InvalidOperationException($"Верхняя частота пропускания fph должна быть ниже верхней частоты среза fsh\r\n  dt={dt}\r\n  fsl={fsl}\r\n  fpl={fpl}\r\n  fph={fph}\r\n  fsh={fsh}")
            {
                Data =
                {
                    { nameof(dt), dt },
                    { "fd", 1/dt },
                    { nameof(fsl), fsl },
                    { nameof(fpl), fpl },
                    { nameof(fph), fph },
                    { nameof(fsh), fsh },
                }
            };

        if (fsh >= 1 / dt / 2)
            throw new InvalidOperationException($"Верхняя частота среза fsh должна быть ниже половины частоты дискретизации fd={1 / dt} (1 / (dt={dt}))\r\n  dt={dt}\r\n  fsl={fsl}\r\n  fpl={fpl}\r\n  fph={fph}\r\n  fsh={fsh}")
            {
                Data =
                {
                    { nameof(dt), dt },
                    { "fd", 1/dt },
                    { nameof(fsl), fsl },
                    { nameof(fpl), fpl },
                    { nameof(fph), fph },
                    { nameof(fsh), fsh },
                }
            };
    }

    private static Specification GetSpecification(
        double dt,
        double fsl,
        double fpl,
        double fph,
        double fsh,
        double Gp,
        double Gs)
    {
        CheckFrequencies(dt, fsl, fpl, fph, fsh);

        if (Gp <= Gs)
            throw new ArgumentOutOfRangeException(
                nameof(Gp), Gp, $"Уровень АЧХ в полосе пропускания Gp={Gp} был меньше, либо равен уровню АЧХ в полосе заграждения Gs={Gs}");

        var Fsl = ToAnalogFrequency(fsl, dt);
        var Fpl = ToAnalogFrequency(fpl, dt);
        var Fph = ToAnalogFrequency(fph, dt);
        var Fsh = ToAnalogFrequency(fsh, dt);

        var Wsl = Consts.pi2 * Fsl;
        var Wpl = Consts.pi2 * Fpl;
        var Wph = Consts.pi2 * Fph;
        var Wsh = Consts.pi2 * Fsh;

        var Wc = Wpl * Wph;
        var dW = Wph - Wpl;

        // Выбор опорной частоты
        // Если   Wc / Wsh > Wsl
        // то есть      Wc > Wsl*Wsh
        // то есть Wsl*Wsh > Wsl*Wsh
        // то есть центральная частота по границам подавления > центральной частоты по границам пропускания
        // то выбираем в качестве опорной частоты выбираем верхнюю границу пропускания
        // иначе, выбираем нижнюю границу пропускания
        var Ws = Wc / Wsh > Wsl
            ? Wsh
            : Wsl;
        var Fp = Abs(dW * Ws / (Wc - Ws.Pow2())) / Consts.pi2;
        const double Fs = 1 / Consts.pi2;

        // Для передачи информации о граничных частотах в спецификацию аналогвого прототипа перечситываем частоты цифрового фильтра обратно
        var fp = ToDigitalFrequency(Fp, dt);
        var fs = ToDigitalFrequency(Fs, dt);

        return new Specification(dt, fp, fs, Gp, Gs);
    }

    private static (double[] A, double[] B) Initialize(double fsl, double fpl, double fph, double fsh, Specification Spec)
    {
        // Пересчитываем аналоговые частоты полосы заграждения в цифровые
        var dt = Spec.dt;

        CheckFrequencies(dt, fsl, fpl, fph, fsh);

        var Wpl = Consts.pi2 * ToAnalogFrequency(fpl, dt);
        var Wph = Consts.pi2 * ToAnalogFrequency(fph, dt);

        var kEps = 1 / Spec.kEps;
        var kW = 1 / Spec.kW;

        var Kw = FullEllipticIntegral(kW);
        var Tw = FullEllipticIntegralComplimentary(kW);
        var K_eps = FullEllipticIntegral(kEps);
        var T_eps = FullEllipticIntegralComplimentary(kEps);

        var N = (int)Ceiling(T_eps * Kw / (K_eps * Tw));
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var (zeros, poles) = GetNormedZerosPoles(N, Spec.EpsP, Spec.EpsS);

        var Fpl = ToAnalogFrequency(fpl, Spec.dt);
        var Fph = ToAnalogFrequency(fph, Spec.dt);

        var ppf_zeros = TransformToBandPassW(zeros, Wpl, Wph);
        var ppf_poles = TransformToBandPassW(poles, Wpl, Wph);

        var z_zeros_enum = ToZ(ppf_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = ToZArray(ppf_poles, dt);

        var Fp0 = (Fpl * Fph).Sqrt();
        var ffp0 = ToDigitalFrequency(Fp0, dt);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);

        var norm_0 = z_zeros.Multiply(z => z0 - z);
        var norm_p = z_poles.Multiply(z => z0 - z);

        var g_norm = N.IsEven() 
            ? Spec.Gp * (norm_p / norm_0).Abs 
            : (z0 * norm_p / norm_0).Abs;

        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b.Re * g_norm);
        var A = GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    public EllipticBandPass(
        double dt,
        double fsl,
        double fpl,
        double fph,
        double fsh,
        double Gp = 0.891250938,
        double Gs = 0.01)
        : this(fsl, fpl, fph, fsh, GetSpecification(dt, fsl, fpl, fph, fsh, Gp, Gs)) { }

    private EllipticBandPass(double fsl, double fpl, double fph, double fsh, Specification Spec)
        : this(Initialize(fsl, fpl, fph, fsh, Spec), Spec) { }

    private EllipticBandPass((double[] A, double[] B) Polynoms, Specification Spec) : this(Polynoms.B, Polynoms.A, Spec) { }

    private EllipticBandPass(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
}