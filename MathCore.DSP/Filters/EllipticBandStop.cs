using System.Data;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters;

/// <summary>Эллиптический полосо-заграждающий фильтр</summary>
public class EllipticBandStop : EllipticFilter
{
    /// <summary>Формирование спецификации фильтра</summary>
    /// <param name="dt">Период дискретизации сигнала</param>
    /// <param name="fpl">Частота нижней границы полосы пропускания</param>
    /// <param name="fsl">Частота нижней границы полосы заграждения</param>
    /// <param name="fsh">Частота верхней границы полосы заграждения</param>
    /// <param name="fph">Частота верхней границы полосы пропускания</param>
    /// <param name="Gp">Уровень АЧХ в полосе пропускания (в единицах 0..1)</param>
    /// <param name="Gs">Неравномерность АЧХ в полосе заграждения (в единицах 0..1)</param>
    /// <returns>Спецификация фильтра</returns>
    private static Specification GetSpecification(
        double dt,
        double fpl,
        double fsl,
        double fsh,
        double fph,
        double Gp,
        double Gs)
    {
        if (Gp <= Gs)
            throw new ArgumentOutOfRangeException(
                nameof(Gp), Gp, $"Уровень АЧХ в полосе пропускания Gp={Gp} был меньше, либо равен уровню АЧХ в полосе заграждения Gs={Gs}");

        var Fpl = ToDigitalFrequency(fpl, dt);
        var Fsl = ToDigitalFrequency(fsl, dt);
        var Fsh = ToDigitalFrequency(fsh, dt);
        var Fph = ToDigitalFrequency(fph, dt);

        var Wpl = Consts.pi2 * Fpl;
        var Wsl = Consts.pi2 * Fsl;
        var Wsh = Consts.pi2 * Fsh;
        var Wph = Consts.pi2 * Fph;

        var Wc = Wsl * Wsh;
        var dW = Wsh - Wsl;

        // Выбор опорной частоты
        // Если   Wc / Wph > Wpl
        // то есть      Wc > Wpl*Wph
        // то есть Wsl*Wsh > Wpl*Wph
        // то есть центральная частота по границам подавления > центральной частоты по границам пропускания
        // то выбираем в качестве опорной частоты выбираем верхнюю границу пропускания
        // иначе, выбираем нижнюю границу пропускания
        var Wp = Wc / Wph > Wpl
            ? Wph
            : Wpl;
        var Fp = Math.Abs(dW * Wp / (Wc - Wp.Pow2())) / Consts.pi2;
        const double Fs = 1 / Consts.pi2;

        // Для передачи информации о граничных частотах в спецификацию аналогвого прототипа перечситываем частоты цифрового фильтра обратно
        var fp = ToAnalogFrequency(Fp, dt);
        var fs = ToAnalogFrequency(Fs, dt);

        return new Specification(dt, fp, fs, Gp, Gs);
    }

    /// <summary>Расчёт коэффициентов полиномов числителя из знаменателя передаточной функции фильтра</summary>
    /// <param name="fsl">Нижняя частота полосы подавления</param>
    /// <param name="fsh">Верхняя частота полосы подавления</param>
    /// <param name="Spec">Спецификация фильтра</param>
    /// <returns>Кортеж, содержащий массивы A - коэффициенты полинома знаменателя и B - коэффициенты полинома числителя</returns>
    private static (double[] A, double[] B) Initialize(double fsl, double fsh, Specification Spec)
    {
        var kW = 1 / Spec.kW;
        var kEps = 1 / Spec.kEps;

        var N = (int)Math.Ceiling(T(kEps) * K(kW) / (K(kEps) * T(kW))); // Порядок фильтра

        // Нули и полюса аналогового прототипа - эллиптического ФНЧ
        var (zeros, poles) = GetNormedZerosPoles(N, Spec.EpsP, Spec.EpsS, Spec.Wp);

        // Пересчитываем аналоговые частоты полосы заграждения в цифровые
        var dt = Spec.dt;
        var Fsl = ToDigitalFrequency(fsl, dt);
        var Fsh = ToDigitalFrequency(fsh, dt);

        // Переносим нули и полюса аналогового нормированного ФНЧ в полосы ПЗФ
        var pzf_zeros = TransformToBandStop(zeros, Fsl, Fsh);
        var pzf_poles = TransformToBandStop(poles, Fsl, Fsh);

        // Если фильтр нечётный, то надо добавить ещё одну пару нулей на центральной частоте полосы заграждения
        if (N.IsOdd())
        {
            var wc_sqrt = Consts.pi2 * (fsl * fsh).Sqrt(); // странно что не Fsl и Fsh
            pzf_zeros = pzf_zeros.AppendLast(
                (0, +wc_sqrt),
                (0, -wc_sqrt));
        }

        // Преобразуем аналоговые нули и полюса в нули и полюса цифрового фильтра с помощью Z-преобразования
        var z_zeros = ToZArray(pzf_zeros, dt);
        var z_poles = ToZArray(pzf_poles, dt);

        // Вычисляем коэффициент нормировки фильтра на нулевой частоте 
        var g_norm = (N.IsOdd() ? 1 : Spec.Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToArray(b => b.Re * g_norm);
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="dt">Период дискретизации цифрового сигнала</param>
    /// <param name="fpl">Нижняя граница полосы пропускания</param>
    /// <param name="fsl">Нижняя граница полосы подавления</param>
    /// <param name="fsh">Верхняя граница полосы подавления</param>
    /// <param name="fph">Верхняя граница полосы пропускания</param>
    /// <param name="Gp">Уровень АЧХ в полосе пропускания (по умолчанию -1дБ)</param>
    /// <param name="Gs">Уровень АЧХ в полосе подавления (по умолчанию -30дБ)</param>
    public EllipticBandStop(
        double dt,
        double fpl,
        double fsl,
        double fsh,
        double fph,
        double Gp = 0.89125093813374556,
        double Gs = 0.01)
        : this(fsl, fsh, GetSpecification(dt, fpl, fsl, fsh, fph, Gp, Gs)) { }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="fsl">Нижняя граница полосы подавления</param>
    /// <param name="fsh">Верхняя граница полосы подавления</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private EllipticBandStop(double fsl, double fsh, Specification Spec) : this(Initialize(fsl, fsh, Spec), Spec) { }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="Polynoms">Кортеж с коэффициентами полиномов знаменателя и числителя функции фильтра</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private EllipticBandStop((double[] A, double[] B) Polynoms, Specification Spec) : this(Polynoms.B, Polynoms.A, Spec) { }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="B">Коэффициенты полинома числителя</param>
    /// <param name="A">Коэффициенты полинома знаменателя</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private EllipticBandStop(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
}