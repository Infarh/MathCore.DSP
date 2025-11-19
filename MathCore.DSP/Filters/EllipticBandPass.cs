using System.Diagnostics;

using MathCore.Extensions;

using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions.EllipticJacobi;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters;

/// <summary>Полосовой эллиптический фильтр</summary>
public class EllipticBandPass : EllipticFilter
{
    /// <summary>Проверяет корректность частот полосы фильтра</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="fsl">Нижняя частота среза</param>
    /// <param name="fpl">Нижняя частота пропускания</param>
    /// <param name="fph">Верхняя частота пропускания</param>
    /// <param name="fsh">Верхняя частота среза</param>
    private static void CheckFrequencies(double dt, double fsl, double fpl, double fph, double fsh)
    {
        if (dt <= 0)
            throw new InvalidOperationException($"Период дискретизации dt={dt} не может быть меньше, либо равен нулю")
                .WithData(nameof(dt), dt)
                .WithData("fd", 1 / dt)
                .WithData(nameof(fsl), fsl)
                .WithData(nameof(fpl), fpl)
                .WithData(nameof(fph), fph)
                .WithData(nameof(fsh), fsh);
        if (1 / dt == 0)
            throw new InvalidOperationException("Частота дискретизации не может быть равна нулю")
                .WithData(nameof(dt), dt)
                .WithData("fd", 1 / dt)
                .WithData(nameof(fsl), fsl)
                .WithData(nameof(fpl), fpl)
                .WithData(nameof(fph), fph)
                .WithData(nameof(fsh), fsh);
        if (fsl >= fpl)
            throw new InvalidOperationException($"""
                Нижняя частота среза fsl должна быть ниже нижней частоты пропускания fpl
                  dt={dt}
                  fsl={fsl}
                  fpl={fpl}
                  fph={fph}
                  fsh={fsh}
                """)
                .WithData(nameof(dt), dt)
                .WithData("fd", 1 / dt)
                .WithData(nameof(fsl), fsl)
                .WithData(nameof(fpl), fpl)
                .WithData(nameof(fph), fph)
                .WithData(nameof(fsh), fsh);
        if (fpl >= fph)
            throw new InvalidOperationException($"""
                Нижняя частота пропускания fpl должна быть ниже верхней частоты пропускания fph
                  dt={dt}
                  fsl={fsl}
                  fpl={fpl}
                  fph={fph}
                  fsh={fsh}
                """)
                .WithData(nameof(dt), dt)
                .WithData("fd", 1 / dt)
                .WithData(nameof(fsl), fsl)
                .WithData(nameof(fpl), fpl)
                .WithData(nameof(fph), fph)
                .WithData(nameof(fsh), fsh);
        if (fph >= fsh)
            throw new InvalidOperationException($"""
                Верхняя частота пропускания fph должна быть ниже верхней частоты среза fsh
                  dt={dt}
                  fsl={fsl}
                  fpl={fpl}
                  fph={fph}
                  fsh={fsh}
                """)
                .WithData(nameof(dt), dt)
                .WithData("fd", 1 / dt)
                .WithData(nameof(fsl), fsl)
                .WithData(nameof(fpl), fpl)
                .WithData(nameof(fph), fph)
                .WithData(nameof(fsh), fsh);
        if (fsh >= 1 / dt / 2)
            throw new InvalidOperationException($"""
                Верхняя частота среза fsh должна быть ниже половины частоты дискретизации fd={1 / dt} (1 / (dt={dt}))
                  dt={dt}
                  fsl={fsl}
                  fpl={fpl}
                  fph={fph}
                  fsh={fsh}
                """)
                .WithData(nameof(dt), dt)
                .WithData("fd", 1 / dt)
                .WithData(nameof(fsl), fsl)
                .WithData(nameof(fpl), fpl)
                .WithData(nameof(fph), fph)
                .WithData(nameof(fsh), fsh);
    }

    /// <summary>Формирует спецификацию фильтра</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="fsl">Нижняя частота среза</param>
    /// <param name="fpl">Нижняя частота пропускания</param>
    /// <param name="fph">Верхняя частота пропускания</param>
    /// <param name="fsh">Верхняя частота среза</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
    /// <returns>Спецификация фильтра</returns>
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

        var Fsl = ToDigitalFrequency(fsl, dt);
        var Fpl = ToDigitalFrequency(fpl, dt);
        var Fph = ToDigitalFrequency(fph, dt);
        var Fsh = ToDigitalFrequency(fsh, dt);

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
        var fp = ToAnalogFrequency(Fp, dt);
        var fs = ToAnalogFrequency(Fs, dt);

        return new(dt, fp, fs, Gp, Gs);
    }

    /// <summary>Вычисляет коэффициенты фильтра</summary>
    /// <param name="fsl">Нижняя частота среза</param>
    /// <param name="fpl">Нижняя частота пропускания</param>
    /// <param name="fph">Верхняя частота пропускания</param>
    /// <param name="fsh">Верхняя частота среза</param>
    /// <param name="Spec">Спецификация фильтра</param>
    /// <returns>Кортеж массивов коэффициентов знаменателя и числителя</returns>
    private static (double[] A, double[] B) Initialize(double fsl, double fpl, double fph, double fsh, Specification Spec)
    {
        // Пересчитываем аналоговые частоты полосы заграждения в цифровые
        var dt = Spec.dt;

        CheckFrequencies(dt, fsl, fpl, fph, fsh);

        var Wpl = Consts.pi2 * ToDigitalFrequency(fpl, dt);
        var Wph = Consts.pi2 * ToDigitalFrequency(fph, dt);

        var kEps = 1 / Spec.kEps;
        var kW = 1 / Spec.kW;

        var Kw = FullEllipticIntegral(kW);
        var Tw = FullEllipticIntegralComplimentary(kW);
        var K_eps = FullEllipticIntegral(kEps);
        var T_eps = FullEllipticIntegralComplimentary(kEps);

        var N = (int)Ceiling(T_eps * Kw / (K_eps * Tw));
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var (zeros, poles) = GetNormedZerosPoles(N, Spec.EpsP, Spec.EpsS);

        var Fpl = ToDigitalFrequency(fpl, Spec.dt);
        var Fph = ToDigitalFrequency(fph, Spec.dt);

        var ppf_zeros = TransformToBandPassW(zeros, Wpl, Wph);
        var ppf_poles = TransformToBandPassW(poles, Wpl, Wph);

        var z_zeros_enum = ToZ(ppf_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendFirst(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = ToZArray(ppf_poles, dt);

        var Fp0 = (Fpl * Fph).Sqrt();
        var ffp0 = ToAnalogFrequency(Fp0, dt);
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

    /// <summary>Создаёт полосовой эллиптический фильтр</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="fsl">Нижняя частота среза</param>
    /// <param name="fpl">Нижняя частота пропускания</param>
    /// <param name="fph">Верхняя частота пропускания</param>
    /// <param name="fsh">Верхняя частота среза</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
    public EllipticBandPass(
        double dt,
        double fsl,
        double fpl,
        double fph,
        double fsh,
        double Gp = 0.891250938,
        double Gs = 0.01)
        : this(fsl, fpl, fph, fsh, GetSpecification(dt, fsl, fpl, fph, fsh, Gp, Gs)) { }

    /// <summary>Создаёт полосовой эллиптический фильтр по спецификации</summary>
    /// <param name="fsl">Нижняя частота среза</param>
    /// <param name="fpl">Нижняя частота пропускания</param>
    /// <param name="fph">Верхняя частота пропускания</param>
    /// <param name="fsh">Верхняя частота среза</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private EllipticBandPass(double fsl, double fpl, double fph, double fsh, Specification Spec)
        : this(Initialize(fsl, fpl, fph, fsh, Spec), Spec) { }

    /// <summary>Создаёт полосовой эллиптический фильтр по массивам коэффициентов</summary>
    /// <param name="Polynoms">Кортеж массивов коэффициентов</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private EllipticBandPass((double[] A, double[] B) Polynoms, Specification Spec) : this(Polynoms.B, Polynoms.A, Spec) { }

    /// <summary>Создаёт полосовой эллиптический фильтр по массивам коэффициентов</summary>
    /// <param name="B">Коэффициенты числителя</param>
    /// <param name="A">Коэффициенты знаменателя</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private EllipticBandPass(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
}