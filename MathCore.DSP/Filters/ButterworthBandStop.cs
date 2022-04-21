using static System.Math;
using static MathCore.Polynom.Array;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters;

public class ButterworthBandStop : ButterworthFilter
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
        double Gp = 0.891250938,
        double Gs = 0.031622777)
    {
        if (Gp <= Gs)
            throw new ArgumentOutOfRangeException(
                nameof(Gp), Gp, $"Уровень АЧХ в полосе пропускания Gp={Gp} был меньше, либо равен уровню АЧХ в полосе заграждения Gs={Gs}");

        var Fpl = ToAnalogFrequency(fpl, dt);
        var Fsl = ToAnalogFrequency(fsl, dt);
        var Fsh = ToAnalogFrequency(fsh, dt);
        var Fph = ToAnalogFrequency(fph, dt);

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
        var W0 = Abs(dW * Wp / (Wc - Wp.Pow2())); // пересчитываем выбранную границу в нижнюю границу пропускания АЧХ аналогового прототипа
        const double W1 = 1;                           // верхняя граница АЧХ аналогового прототипа будет всегда равна 1 рад/с
        var Fp = W0 / Consts.pi2;
        const double Fs = 1 / Consts.pi2;

        // Для передачи информации о граничных частотах в спецификацию аналогвого прототипа перечситываем частоты цифрового фильтра обратно
        var fp = ToDigitalFrequency(Fp, dt);
        var fs = ToDigitalFrequency(Fs, dt);

        return new Specification(dt, fp, fs, Gp, Gs);
    }

    /// <summary>Расчёт коэффициентов полиномов числителя из знаменателя передаточной функции фильтра</summary>
    /// <param name="fpl">Нижняя частота полосы пропускания</param>
    /// <param name="fsl">Нижняя частота полосы подавления</param>
    /// <param name="fsh">Верхняя частота полосы подавления</param>
    /// <param name="fph">Верхняя частота полосы пропускания</param>
    /// <param name="Spec">Спецификация фильтра</param>
    /// <returns>Кортеж, содержащий массивы A - коэффициенты полинома знаменателя и B - коэффициенты полинома числителя</returns>
    private static (double[] A, double[] B) Initialize(double fpl, double fsl, double fsh, double fph, Specification Spec)
    {
        // Пересчитываем аналоговые частоты полосы заграждения в цифровые
        var dt = Spec.dt;

        var Wpl = Consts.pi2 * ToAnalogFrequency(fpl, dt);
        var Wsl = Consts.pi2 * ToAnalogFrequency(fsl, dt);
        var Wsh = Consts.pi2 * ToAnalogFrequency(fsh, dt);
        var Wph = Consts.pi2 * ToAnalogFrequency(fph, dt);

        var Wc = Wsl * Wsh;
        var dW = Wsh - Wsl;
        var sqrtWc = Wc.Sqrt();
        var Wp = Wc / Wph > Wpl
            ? Wph
            : Wpl;
        var W0 = Abs(dW * Wp / (Wc - Wp.Pow2()));

        var N = (int)Ceiling(Log(Spec.kEps) / Log(Spec.kW));
        var poles = GetNormPoles(N, Spec.EpsP, W0);

        // Переносим нули и полюса аналогового нормированного ФНЧ в полосы ПЗФ
        var pzf_zeros = Enumerable.Range(0, 2 * N).Select(i => i % 2 == 0 ? Complex.ImValue(sqrtWc) : Complex.ImValue(-sqrtWc));
        var pzf_poles = TransformToBandStopW(poles, Wsl, Wsh);

        // Преобразуем аналоговые нули и полюса в нули и полюса цифрового фильтра с помощью Z-преобразования
        var z_zeros = ToZArray(pzf_zeros, dt);
        var z_poles = ToZArray(pzf_poles, dt);

        // Вычисляем коэффициент нормировки фильтра на нулевой частоте 
        var G_norm = (N.IsOdd() ? 1 : Spec.Gp)
            / (z_zeros.Multiply(z => 1 - z) / z_poles.Multiply(z => 1 - z)).Abs;

        // Определяем массивы нулей коэффициентов полиномов знаменателя и числителя
        var B = GetCoefficientsInverted(z_zeros).ToArray(b => b.Re * G_norm);
        var A = GetCoefficientsInverted(z_poles).ToRe();

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
    public ButterworthBandStop(
        double dt,
        double fpl,
        double fsl,
        double fsh,
        double fph,
        double Gp = 0.89125093813374556,
        double Gs = 0.031622777) 
        : this(fpl, fsl, fsh, fph, GetSpecification(dt, fpl, fsl, fsh, fph, Gp, Gs)) { }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="fsl">Нижняя граница полосы подавления</param>
    /// <param name="fsh">Верхняя граница полосы подавления</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private ButterworthBandStop(double fpl,double fsl, double fsh, double fph, Specification Spec) 
        : this(Initialize(fpl, fsl, fsh, fph, Spec), Spec) { }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="Polynoms">Кортеж с коэффициентами полиномов знаменателя и числителя функции фильтра</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private ButterworthBandStop((double[] A, double[] B) Polynoms, Specification Spec) : this(Polynoms.B, Polynoms.A, Spec) { }

    /// <summary>Инициализация нового эллиптического полосозаграждающего фильтра (ПЗФ)</summary>
    /// <param name="B">Коэффициенты полинома числителя</param>
    /// <param name="A">Коэффициенты полинома знаменателя</param>
    /// <param name="Spec">Спецификация фильтра</param>
    private ButterworthBandStop(double[] B, double[] A, Specification Spec) : base(B, A, Spec) { }
}