﻿using System.ComponentModel;
using System.Diagnostics;

using static System.Math;

namespace MathCore.DSP.Filters;

public class ChebyshevBandPass : ChebyshevFilter
{
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
        double Gp,
        double Gs)
    {
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
        // то есть Wpl*Wph > Wsl*Wsh
        // то есть центральная частота по границам подавления > центральной частоты по границам пропускания
        // то выбираем в качестве опорной частоты выбираем верхнюю границу пропускания
        // иначе, выбираем нижнюю границу пропускания
        var Ws = Wc / Wsh > Wsl
            ? Wsh
            : Wsl;
        //const double W0 = 1;                           // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var W1 = Abs((Wc - Ws.Pow2()) / (dW * Ws)); // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        const double Fp = 1 / Consts.pi2;
        var Fs = W1 / Consts.pi2;

        // Для передачи информации о граничных частотах в спецификацию аналогвого прототипа перечситываем частоты цифрового фильтра обратно
        var fp = ToAnalogFrequency(Fp, dt);
        var fs = ToAnalogFrequency(Fs, dt);

        return new(dt, fp, fs, Gp, Gs);
    }

    /// <summary>Расчёт коэффициентов полиномов числителя из знаменателя передаточной функции фильтра</summary>
    /// <param name="fpl">Нижняя частота полосы пропускания</param>
    /// <param name="fph">Верхняя частота полосы пропускания</param>
    /// <param name="Spec">Спецификация фильтра</param>
    /// <returns>Кортеж, содержащий массивы A - коэффициенты полинома знаменателя и B - коэффициенты полинома числителя</returns>
    private static (double[] A, double[] B) InitializeI(double fpl, double fph, Specification Spec)
    {
        // Пересчитываем аналоговые частоты полосы заграждения в цифровые
        var dt = Spec.dt;

        var Wpl = Consts.pi2 * ToDigitalFrequency(fpl, dt);
        var Wph = Consts.pi2 * ToDigitalFrequency(fph, dt);

        var N = (int)Ceiling(arcch(Spec.kEps) / arcch(Spec.kW));
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var poles = GetNormedPolesI(N, Spec.EpsP);

        // Переносим нули и полюса аналогового нормированного ФНЧ в полосы ПЗФ
        var pzf_poles = TransformToBandPassW(poles, Wpl, Wph);

        // Преобразуем аналоговые нули и полюса в нули и полюса цифрового фильтра с помощью Z-преобразования
        var z_zeros = Enumerable
           .Repeat(Complex.Real, N)
           .Concat(Enumerable.Repeat(-Complex.Real, N))
           .ToArray();
        var z_poles = ToZArray(pzf_poles, dt);

        var Fpl = Wpl / Consts.pi2;
        var Fph = Wph / Consts.pi2;
        var ffp0 = ToAnalogFrequency((Fpl * Fph).Sqrt(), dt);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);

        // Вычисляем коэффициент нормировки фильтра на нулевой частоте 
        var norm_0 = ((z0 - 1) * (z0 + 1)).Pow(N).Abs;
        var norm_p = z_poles.Multiply(z => z0 - z);

        var g_norm = N.IsEven() 
            ? Spec.Gp * (norm_p / norm_0).Abs 
            : (z0 * norm_p / norm_0).Abs;

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToArray(b => b.Re * g_norm);
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    private static (double[] A, double[] B) InitializeII(double fpl, double fph, Specification Spec)
    {
        // Пересчитываем аналоговые частоты полосы заграждения в цифровые
        var dt = Spec.dt;

        var Wpl = Consts.pi2 * ToDigitalFrequency(fpl, dt);
        var Wph = Consts.pi2 * ToDigitalFrequency(fph, dt);

        var N = (int)Ceiling(arcch(Spec.kEps) / arcch(Spec.kW));
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var (zeros, poles) = GetNormedPolesII(N, Spec.EpsS);

        // Переносим нули и полюса аналогового нормированного ФНЧ в полосы ПЗФ
        var ppf_zeros_enum = TransformToBandPassW(zeros, Wpl, Wph);
        if (N.IsOdd())
            ppf_zeros_enum = ppf_zeros_enum.AppendLast(0);
        var ppf_zeros = ppf_zeros_enum.ToArray();
        var ppf_poles = TransformToBandPassW(poles, Wpl, Wph);

        // Преобразуем аналоговые нули и полюса в нули и полюса цифрового фильтра с помощью Z-преобразования
        var z_zeros_enum = ToZ(ppf_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendLast(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = ToZArray(ppf_poles, dt);

        // Вычисляем коэффициент нормировки фильтра на нулевой частоте 
        var ffp0 = ToAnalogFrequency((Wpl * Wph).Sqrt() / Consts.pi2, dt);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);

        var norm_0 = z_zeros.Multiply(z => z0 - z);
        var norm_p = z_poles.Multiply(z => z0 - z);

        var g_norm = N.IsEven() 
            ? Spec.Gp * (norm_p / norm_0).Re 
            : (norm_p / norm_0).Re;

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToArray(b => b.Re * g_norm);
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    private static (double[] A, double[] B) InitializeIICorrected(double fpl, double fph, Specification Spec)
    {
        // Пересчитываем аналоговые частоты полосы заграждения в цифровые
        var dt = Spec.dt;

        var Wpl = Consts.pi2 * ToDigitalFrequency(fpl, dt);
        var Wph = Consts.pi2 * ToDigitalFrequency(fph, dt);

        var N = (int)Ceiling(arcch(Spec.kEps) / arcch(Spec.kW));
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var (zeros, poles) = GetNormedPolesII(N, Spec.EpsS, Spec.kW);

        // Переносим нули и полюса аналогового нормированного ФНЧ в полосы ПЗФ
        var ppf_zeros_enum = AnalogBasedFilter.TransformToBandPassW(zeros, Wpl, Wph);
        if (N.IsOdd())
            ppf_zeros_enum = ppf_zeros_enum.AppendLast(0);
        var ppf_zeros = ppf_zeros_enum.ToArray();
        var ppf_poles = TransformToBandPassW(poles, Wpl, Wph).ToArray();

        // Преобразуем аналоговые нули и полюса в нули и полюса цифрового фильтра с помощью Z-преобразования
        var z_zeros_enum = ToZ(ppf_zeros, dt);
        if (N.IsOdd())
            z_zeros_enum = z_zeros_enum.AppendLast(-1);
        var z_zeros = z_zeros_enum.ToArray();
        var z_poles = ToZArray(ppf_poles, dt);

        // Вычисляем коэффициент нормировки фильтра на нулевой частоте 
        var ffp0 = ToAnalogFrequency((Wpl * Wph).Sqrt() / Consts.pi2, dt);
        var z0 = Complex.Exp(Consts.pi2 * ffp0 * dt);

        var norm_0 = z_zeros.Multiply(z => z0 - z);
        var norm_p = z_poles.Multiply(z => z0 - z);

        var g_norm = N.IsEven() 
            ? Spec.Gp * (norm_p / norm_0).Abs
            : (z0 * norm_p / norm_0).Abs;

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToArray(b => b.Re * g_norm);
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="dt">Период дискретизации цифрового сигнала</param>
    /// <param name="fsl">Нижняя граница полосы подавления</param>
    /// <param name="fpl">Нижняя граница полосы пропускания</param>
    /// <param name="fph">Верхняя граница полосы пропускания</param>
    /// <param name="fsh">Верхняя граница полосы подавления</param>
    /// <param name="Gp">Уровень АЧХ в полосе пропускания (по умолчанию -1дБ)</param>
    /// <param name="Gs">Уровень АЧХ в полосе подавления (по умолчанию -40дБ)</param>
    /// <param name="Type">Тип фильтра I или II</param>
    public ChebyshevBandPass(
        double dt,
        double fsl,
        double fpl,
        double fph,
        double fsh,
        double Gp = 0.89125093813374556,
        double Gs = 0.01, 
        ChebyshevType Type = ChebyshevType.I)
        : this(fpl, fph, GetSpecification(dt, fsl, fpl, fph, fsh, Gp, Gs), Type) { }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="fpl">Нижняя граница полосы пропускания</param>
    /// <param name="fph">Верхняя граница полосы пропускания</param>
    /// <param name="Spec">Спецификация фильтра</param>
    /// <param name="Type">Тип фильтра I или II</param>
    private ChebyshevBandPass(double fpl, double fph, Specification Spec, ChebyshevType Type)
        : this(Type switch
        {
            ChebyshevType.I => InitializeI(fpl, fph, Spec),
            ChebyshevType.II => InitializeII(fpl, fph, Spec),
            ChebyshevType.IICorrected => InitializeIICorrected(fpl, fph, Spec),
            _ => throw new InvalidEnumArgumentException(nameof(Type), (int)Type, typeof(ChebyshevType))
        }, Spec, Type) { }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="Polynoms">Кортеж с коэффициентами полиномов знаменателя и числителя функции фильтра</param>
    /// <param name="Spec">Спецификация фильтра</param>
    /// <param name="Type">Тип фильтра I или II</param>
    private ChebyshevBandPass((double[] A, double[] B) Polynoms, Specification Spec, ChebyshevType Type) 
        : this(Polynoms.B, Polynoms.A, Spec, Type) { }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="B">Коэффициенты полинома числителя</param>
    /// <param name="A">Коэффициенты полинома знаменателя</param>
    /// <param name="Spec">Спецификация фильтра</param>
    /// <param name="Type">Тип фильтра I или II</param>
    private ChebyshevBandPass(double[] B, double[] A, Specification Spec, ChebyshevType Type) : base(B, A, Spec, Type) { }
}