using System.Diagnostics;

using static System.Math;

using static MathCore.Consts;
using static MathCore.Polynom.Array;

// ReSharper disable InconsistentNaming
// ReSharper disable MemberCanBePrivate.Global
// ReSharper disable UnusedAutoPropertyAccessor.Global

namespace MathCore.DSP.Filters;

public class ButterworthBandPass : ButterworthFilter
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

    /// <summary>Формирование спецификации фильтра</summary>
    /// <param name="dt">Период дискретизации сигнала</param>
    /// <param name="fsl">Частота нижней границы полосы заграждения</param>
    /// <param name="fpl">Частота нижней границы полосы пропускания</param>
    /// <param name="fph">Частота верхней границы полосы пропускания</param>
    /// <param name="fsh">Частота верхней границы полосы заграждения</param>
    /// <param name="Gp">Уровень АЧХ в полосе пропускания (в единицах 0..1)</param>
    /// <param name="Gs">Неравномерность АЧХ в полосе заграждения (в единицах 0..1)</param>
    /// <returns>Спецификация фильтра</returns>
    private static Specification GetSpecification(
        double dt,
        double fsl,
        double fpl,
        double fph,
        double fsh,
        double Gp = 0.891250938,
        double Gs = 0.031622777)
    {
        CheckFrequencies(dt, fsl, fpl, fph, fsh);

        if (Gp <= Gs)
            throw new ArgumentOutOfRangeException(
                nameof(Gp), Gp, $"Уровень АЧХ в полосе пропускания Gp={Gp} был меньше, либо равен уровню АЧХ в полосе заграждения Gs={Gs}");

        var Fsl = ToDigitalFrequency(fsl, dt);
        var Fpl = ToDigitalFrequency(fpl, dt);
        var Fph = ToDigitalFrequency(fph, dt);
        var Fsh = ToDigitalFrequency(fsh, dt);

        var Wsl = pi2 * Fsl;
        var Wpl = pi2 * Fpl;
        var Wph = pi2 * Fph;
        var Wsh = pi2 * Fsh;

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
        var Fp = Abs(dW * Ws / (Wc - Ws.Pow2())) / pi2;
        const double Fs = 1 / pi2;

        // Для передачи информации о граничных частотах в спецификацию аналогвого прототипа перечситываем частоты цифрового фильтра обратно
        var fp = ToAnalogFrequency(Fp, dt);
        var fs = ToAnalogFrequency(Fs, dt);

        return new Specification(dt, fp, fs, Gp, Gs);
    }

    private static (double[] A, double[] B) Initialize(double fpl, double fph, Specification Spec)
    {
        // Пересчитываем аналоговые частоты полосы заграждения в цифровые
        var dt = Spec.dt;

        var Wpl = pi2 * ToDigitalFrequency(fpl, dt);
        var Wph = pi2 * ToDigitalFrequency(fph, dt);

        var dW = Wph - Wpl;

        var N = (int)Ceiling(Log(Spec.kEps) / Log(Spec.kW));
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");
        var poles = GetNormPoles(N, Spec.EpsP);

        var ppf_poles = TransformToBandPassW(poles, Wpl, Wph).ToArray();

        var ppf_zeros_z = Enumerable
           .Repeat(Complex.ReValue(1), N)
           .AppendLast(Enumerable.Repeat(Complex.ReValue(-1), N))
           .ToArray();
        var ppf_poles_z = ToZArray(ppf_poles, Spec.dt);

        var g_norm = Complex.Real;
        var fd2 = 2 * Spec.fd;
        for (var i = 0; i < N; i++)
        {
            var (p1, p2) = (ppf_poles[2 * i], ppf_poles[2 * i + 1]);
            var p = (fd2 - p1) * (fd2 - p2);
            g_norm *= fd2 * dW / p;
        }

        var G_norm = g_norm / Spec.EpsP;

        var B = GetCoefficientsInverted(ppf_zeros_z).ToArray(v => (v * G_norm).Re);
        var A = GetCoefficientsInverted(ppf_poles_z).ToRe();

        return (A, B);
    }

    public ButterworthBandPass(
        double dt,
        double fsl,
        double fpl,
        double fph,
        double fsh,
        double Gp = 0.891250938,
        double Gs = 0.01)
        : this(fsl, fpl, fph, fsh, GetSpecification(dt, fsl, fpl, fph, fsh, Gp, Gs))
    { }

    private ButterworthBandPass(double fsl, double fpl, double fph, double fsh, Specification Spec)
        : this(Initialize(fpl, fph, Spec), Spec) { }

    private ButterworthBandPass((double[] A, double[] B) Polynoms, Specification Spec) : this(Polynoms.B, Polynoms.A, Spec) { }

    private ButterworthBandPass(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
}