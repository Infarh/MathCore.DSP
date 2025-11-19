using System.Runtime.Serialization;

using static MathCore.Consts;
// ReSharper disable InconsistentNaming
// ReSharper disable ParameterHidesMember

namespace MathCore.DSP.Filters;

/// <summary>Цифровой фильтр на основе аналогового прототипа</summary>
//[KnownType(typeof(BesselFilter))]     // Наихудшая аппроксимация прямоугольной АЧХ (близка по форме к гауссовой кривой), но наилучшая форма переходной хар-ки
[KnownType(typeof(ButterworthFilter))]  // Промежуточные качества по прямоугольности АЧХ и переходной хар-ке
[KnownType(typeof(ChebyshevFilter))]    // АЧХ наиболее приближена к прямоугольной - наибольшие пульсации в переходной хор-ке
[KnownType(typeof(EllipticFilter))]     // Максимальная крутизна АЧХ
public abstract class AnalogBasedFilter : IIR
{
    //https://ru.dsplib.org/content/filter_low2low/filter_low2low.html
    /// <summary>Класс преобразований нулей и полюсов аналоговых прототипов</summary>
    public static class Transform
    {
        /// <summary>Преобразование нулей/полюсов фильтра ФНЧ-ФНЧ</summary>
        /// <param name="Z">Нули/полюса прототипа (нормированного ФНЧ)</param>
        /// <param name="wp">Новая частота среза</param>
        /// <returns>Нули/полюса нового фильтра ФНЧ</returns>
        public static IEnumerable<Complex> ToLow(IEnumerable<Complex> Z, double wp) => Z.Select(z => z / wp);

        /// <summary>Преобразование нулей/полюсов фильтра ФНЧ-ФВЧ</summary>
        /// <param name="Z">Нули/полюса прототипа (нормированного ФНЧ)</param>
        /// <param name="wp">Новая частота среза</param>
        /// <returns>Нули/полюса нового фильтра ФВЧ</returns>
        public static IEnumerable<Complex> ToHigh(IEnumerable<Complex> Z, double wp) => Z.Select(z => z * wp);

        /// <summary>Установить значения новых полюсов/нулей</summary>
        /// <param name="p">Исходное значение</param>
        /// <param name="D">Смещение</param>
        /// <param name="p1">Первое значение</param>
        /// <param name="p2">Второе значение</param>
        private static void Set(in Complex p, in Complex D, out Complex p1, out Complex p2) => (p1, p2) = (p + D, p - D);

        /// <summary>Преобразование нулей/полюсов фильтра ФНЧ-ППФ</summary>
        /// <param name="Poles">Полюса прототипа (нормированного ФНЧ)</param>
        /// <param name="Zeros">Нули прототипа (нормированного ФНЧ)</param>
        /// <param name="fpl">Нижняя частота среза ППФ</param>
        /// <param name="fph">Верхняя частота среза ППФ</param>
        /// <returns>Нули/полюса ППФ</returns>
        /// <exception cref="ArgumentException">Если число полюсов == 0</exception>
        /// <exception cref="ArgumentException">Если число нулей больше числа полюсов</exception>
        public static (Complex[] Poles, Complex[] Zeros) ToBandPass(
            Complex[] Poles,
            Complex[] Zeros,
            double fpl, double fph)
        {
            if (Poles is null) throw new ArgumentNullException(nameof(Poles));
            if (Zeros is null) throw new ArgumentNullException(nameof(Zeros));
            if (Poles.Length == 0) throw new ArgumentException("Размер вектора полюсов должен быть больше 0", nameof(Poles));
            if (Zeros.Length > Poles.Length) throw new ArgumentException("Число нулей не должна быть больше числа полюсов", nameof(Zeros));

            var (wpl, wph) = (pi2 * fpl, pi2 * fph);
            var (dw05, wc2) = ((wph - wph) / 2, wpl * wph);

            // На каждый исходный полюс формируется пара новых полюсов
            // Число нулей равно удвоенному числу исходных нулей
            //      + ноль в 0 кратности, равной разности числа полюсов и нулей
            var poles = new Complex[Poles.Length * 2];
            var zeros = new Complex[Zeros.Length * 2 + (Poles.Length - Zeros.Length)];


            for (var i = 0; i < Poles.Length; i++)
            {
                var pdw = dw05 * Poles[i];
                Set(pdw, Complex.Sqrt(pdw.Power - wc2),
                    out poles[2 * i],
                    out poles[2 * i + 1]);
            }

            for (var i = 0; i < Zeros.Length; i++)
            {
                var pdw = dw05 * Zeros[i];
                Set(pdw, Complex.Sqrt(pdw.Power - wc2),
                    out zeros[2 * i],
                    out zeros[2 * i + 1]);
            }

            return (poles, zeros);
        }

        /// <summary>Преобразование нулей/полюсов фильтра ФНЧ-ПЗФ</summary>
        /// <param name="Poles">Полюса прототипа (нормированного ФНЧ)</param>
        /// <param name="Zeros">Нули прототипа (нормированного ФНЧ)</param>
        /// <param name="fpl">Нижняя частота среза ПЗФ</param>
        /// <param name="fph">Верхняя частота среза ПЗФ</param>
        /// <returns>Нули/полюса ППФ</returns>
        /// <exception cref="ArgumentException">Если число полюсов == 0</exception>
        /// <exception cref="ArgumentException">Если число нулей больше числа полюсов</exception>
        public static (Complex[] Poles, Complex[] Zeros) ToBandStop(
            Complex[] Poles,
            Complex[] Zeros,
            double fpl, double fph)
        {
            if (Poles is null) throw new ArgumentNullException(nameof(Poles));
            if (Zeros is null) throw new ArgumentNullException(nameof(Zeros));
            var count_p = Poles.Length;
            var count_0 = Zeros.Length;
            if (count_p == 0) throw new ArgumentException("Размер вектора полюсов должен быть больше 0", nameof(Poles));
            if (count_0 > count_p) throw new ArgumentException("Число нулей не должна быть больше числа полюсов", nameof(Zeros));

            var (wpl, wph) = (pi2 * fpl, pi2 * fph);
            var (dw05, wc2) = ((wph - wph) / 2, wpl * wph);

            // На каждый исходный полюс формируется пара новых полюсов
            // Число нулей равно удвоенному числу исходных нулей
            //      + ноль в 0 кратности, равной разности числа полюсов и нулей
            var poles = new Complex[count_p * 2];
            var zeros = new Complex[poles.Length];

            //static void Set(in Complex p, in Complex D, out Complex p1, out Complex p2)
            //{
            //    p1 = p + D;
            //    p2 = p - D;
            //}

            for (var i = 0; i < count_p; i++)
            {
                var pdw = dw05 / Poles[i];
                Set(pdw, Complex.Sqrt(pdw.Power - wc2),
                    out poles[2 * i],
                    out poles[2 * i + 1]);
            }

            for (var i = 0; i < count_0; i++)
            {
                var pdw = dw05 / Zeros[i];
                Set(pdw, Complex.Sqrt(pdw.Power - wc2),
                    out zeros[2 * i],
                    out zeros[2 * i + 1]);
            }

            var wc = wc2.Sqrt();
            for (var i = 0; i < count_p - count_0; i++)
            {
                zeros[2 * (count_0 + i)] = new(0, wc);
                zeros[2 * (count_0 + i) + 1] = new(0, -wc);
            }

            return (poles, zeros);
        }
    }

    /// <summary>Спецификация фильтра</summary>
    public readonly ref struct Specification
    {
        /// <summary>Период дискретизации</summary>
        public double dt { get; }

        /// <summary>Частота дискретизации</summary>
        public double fd => 1 / dt;

        /// <summary>Граничная частота полосы пропускания</summary>
        public double fp { get; }
        /// <summary>Граничная частота полосы заграждения</summary>
        public double fs { get; }

        /// <summary>Граничная циклическая частота полосы пропускания</summary>
        public double wp { get; }
        /// <summary>Граничная циклическая частота полосы заграждения</summary>
        public double ws { get; }

        /// <summary>Граничная частота полосы пропускания цифрового фильтра</summary>
        public double Fp { get; }
        /// <summary>Граничная частота полосы заграждения цифрового фильтра</summary>
        public double Fs { get; }

        /// <summary>Граничная циклическая частота полосы пропускания цифрового фильтра</summary>
        public double Wp { get; }
        /// <summary>Граничная циклическая частота полосы заграждения цифрового фильтра</summary>
        public double Ws { get; }

        /// <summary>Коэффициент передачи в полосе пропускания</summary>
        public double Gp { get; }
        /// <summary>Коэффициент передачи в полосе подавления</summary>
        public double Gs { get; }

        /// <summary>Коэффициент подавления в полосе пропускания (дБ)</summary>
        public double Rp { get; }
        /// <summary>Коэффициент подавления в полосе заграждения (дБ)</summary>
        public double Rs { get; }

        /// <summary>Неоднородность АЧХ в полосе пропускания</summary>
        public double EpsP { get; }
        /// <summary>Неоднородность АЧХ в полосе заграждения</summary>
        public double EpsS { get; }

        /// <summary>Отношение неоднородностей коэффициентов передачи заграждения и пропускания</summary>
        public double kEps => EpsS / EpsP;
        /// <summary>Отношение частотных полос заграждения и пропускания</summary>
        public double kw => fs / fp;
        /// <summary>Отношение частотных полос заграждения и пропускания цифрового фильтра</summary>
        public double kW => Ws / Wp;

        /// <summary>Инициализация новой спецификации фильтра</summary>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="fp">Граничная частота полосы пропускания</param>
        /// <param name="fs">Граничная частота полосы заграждения</param>
        /// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
        /// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
        /// <exception cref="InvalidOperationException">Если частота пропускания больше, либо равна частоте подавления</exception>
        /// <exception cref="InvalidOperationException">Если частота пропускания больше, либо равна половине частоты дискретизации</exception>
        public Specification(double dt, double fp, double fs, double Gp, double Gs)
        {
            if (!(fp < 1 / (2 * dt))) throw new InvalidOperationException("Частота пропускания должен быть меньше половины полосы дискретизации");
            if (fp < 0) throw new ArgumentOutOfRangeException(nameof(fp), fp, "Частота среза не может быть меньше 0");
            if (fs < 0) throw new ArgumentOutOfRangeException(nameof(fs), fs, "Частота полосы подавления не может быть меньше 0");
            if (Gp < 0) throw new ArgumentOutOfRangeException(nameof(Gp), Gp, "Значение уровня АЧХ в полосе пропускания не может быть отрицательной величиной");
            if (Gs < 0) throw new ArgumentOutOfRangeException(nameof(Gs), Gs, "Значение уровня АЧХ в полосе заграждения не может быть отрицательной величиной");
            if (Gp > 1) throw new ArgumentOutOfRangeException(nameof(Gp), Gp, "Значение уровня АЧХ в полосе пропускания не может быть больше 1");
            if (Gs > 1) throw new ArgumentOutOfRangeException(nameof(Gs), Gs, "Значение уровня АЧХ в полосе заграждения не может быть больше 1");
            if (Gp <= Gs) throw new InvalidOperationException($"Уровень АЧХ в полосе пропускания Gp должен быть больше уровня в полосе заграждения Gs\r\n  Gp={Gp}\r\n  Gs={Gs}");

            this.dt = dt;
            this.fp = fp;
            this.fs = fs;
            this.Gp = Gp;
            this.Gs = Gs;

            Rp = -Gp.In_dB();
            Rs = -Gs.In_dB();

            EpsP = (1 / (Gp * Gp) - 1).Sqrt();
            //EpsP = (10d.Pow(Rp / 10) - 1).Sqrt();
            EpsS = (10d.Pow(Rs / 10) - 1).Sqrt();

            //var tEpsP = (1 / Gp.Pow2() - 1).Sqrt();
            //var tEpsS = (1 / Gs.Pow2() - 1).Sqrt();

            //var tEpsP = (1 - Gp*Gp).Sqrt() / (Gp*Gp);
            //var tEpsS = (1 - Gs*Gs).Sqrt() / (Gs*Gs);

            Fp = ToDigitalFrequency(fp, dt);  // Частота пропускания аналогового прототипа
            Fs = ToDigitalFrequency(fs, dt);  // Частота подавления аналогового пропита

            (wp, ws) = (pi2 * fp, pi2 * fs);
            (Wp, Ws) = (pi2 * Fp, pi2 * Fs);
        }

        /// <summary>Деконструктор спецификации фильтра</summary>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="fp">Граничная частота полосы пропускания</param>
        /// <param name="fs">Граничная частота полосы заграждения</param>
        /// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
        /// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
        public void Deconstruct(out double dt, out double fp, out double fs, out double Gp, out double Gs) =>
            (dt, fp, fs, Gp, Gs) = (this.dt, this.fp, this.fs, this.Gp, this.Gs);
    }

    /// <summary>Создать спецификацию фильтра</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="fp">Граничная частота полосы пропускания</param>
    /// <param name="fs">Граничная частота полосы заграждения</param>
    /// <param name="Gp">Коэффициент передачи в полосе пропускания</param>
    /// <param name="Gs">Коэффициент передачи в полосе заграждения</param>
    /// <returns>Спецификация фильтра</returns>
    public static Specification GetSpecification(double dt, double fp, double fs, double Gp, double Gs) => new(dt, fp, fs, Gp, Gs);

    /// <summary>Период дискретизации</summary>
    public double dt { get; }

    /// <summary>Частота дискретизации</summary>
    public double fd => 1 / dt;

    /// <summary>Граничная частота полосы пропускания</summary>
    public double fp { get; }

    /// <summary>Граничная частота полосы подавления</summary>
    public double fs { get; }

    /// <summary>Затухание в полосе пропускания</summary>
    public double Gp { get; }

    /// <summary>Затухание в полосе подавления</summary>
    public double Gs { get; }

    /// <summary>Спецификация фильтра</summary>
    public Specification Spec => new(dt, fp, fs, Gp, Gs);

    /// <summary>Инициализация параметров цифрового фильтра на базе аналогового прототипа</summary>
    /// <param name="B">Коэффициенты полинома числителя</param>
    /// <param name="A">Коэффициенты полинома знаменателя</param>
    /// <param name="Spec">Спецификация фильтра</param>
    protected AnalogBasedFilter(double[] B, double[] A, Specification Spec) : base(B, A) => (dt, fp, fs, Gp, Gs) = Spec;

    /// <inheritdoc/>
    public override Complex FrequencyResponse(double f) => base.FrequencyResponse(f / fd);

    /// <inheritdoc/>
    public override Complex FrequencyResponse(double f, double dt) => base.FrequencyResponse(f * dt);

    /// <summary>Метод преобразования нулей и полюсов нормированного ФНЧ в нули и полюса ППФ</summary>
    /// <param name="Normed">Нормированные нули и полюса ФНЧ</param>
    /// <param name="fmin">Нижняя частота среза</param>
    /// <param name="fmax">Верхняя частота среза</param>
    /// <returns>Нули и полюса ППФ</returns>
    public static IEnumerable<Complex> TransformToBandPass(IEnumerable<Complex> Normed, double fmin, double fmax) =>
        TransformToBandPassW(Normed, pi2 * fmin, pi2 * fmax);

    /// <summary>Метод преобразования нулей и полюсов нормированного ФНЧ в нули и полюса ППФ</summary>
    /// <param name="Normed">Нормированные нули и полюса ФНЧ</param>
    /// <param name="w_min">Нижняя частота среза</param>
    /// <param name="w_max">Верхняя частота среза</param>
    /// <returns>Нули и полюса ППФ</returns>
    public static IEnumerable<Complex> TransformToBandPassW(IEnumerable<Complex> Normed, double w_min, double w_max)
    {
        var (dw05, wc2) = ((w_max - w_min) / 2, w_min * w_max);

        foreach (var p in Normed)
        {
            var pdw = dw05 * p;
            var sqrt = Complex.Sqrt(pdw.Pow2() - wc2);
            yield return pdw + sqrt;
            yield return pdw - sqrt;
        }
    }

    /// <summary>Метод преобразования нулей и полюсов нормированного ФНЧ в нули и полюса ПЗФ</summary>
    /// <param name="Normed">Нормированные нули и полюса ФНЧ</param>
    /// <param name="fmin">Нижняя частота среза</param>
    /// <param name="fmax">Верхняя частота среза</param>
    /// <returns>Нули и полюса ПЗФ</returns>
    public static IEnumerable<Complex> TransformToBandStop(IEnumerable<Complex> Normed, double fmin, double fmax) =>
        TransformToBandStopW(Normed, pi2 * fmin, pi2 * fmax);

    /// <summary>Метод преобразования нулей и полюсов нормированного ФНЧ в нули и полюса ПЗФ</summary>
    /// <param name="Normed">Нормированные нули и полюса ФНЧ</param>
    /// <param name="w_min">Нижняя частота среза</param>
    /// <param name="w_max">Верхняя частота среза</param>
    /// <returns>Нули и полюса ПЗФ</returns>
    public static IEnumerable<Complex> TransformToBandStopW(IEnumerable<Complex> Normed, double w_min, double w_max)
    {
        var (dw05, wc2) = ((w_max - w_min) / 2, w_min * w_max);

        foreach (var p in Normed)
        {
            var pdw = dw05 / p;
            var sqrt = Complex.Sqrt(pdw.Pow2() - wc2);
            yield return pdw + sqrt;
            yield return pdw - sqrt;
        }
    }

    /// <summary>Преобразование нулей/полюсов нормированного ФНЧ в нули/полюса ФНЧ с другой частотой среза</summary>
    /// <param name="Normed">Нормированные нули/полюса ФНЧ</param>
    /// <param name="fp">Новая частота среза</param>
    /// <returns>Нули/полюса нового ФНЧ</returns>
    public static IEnumerable<Complex> TransformToLowPass(IEnumerable<Complex> Normed, double fp) =>
        TransformToLowPassW(Normed, pi2 * fp);

    /// <summary>Преобразование нулей/полюсов нормированного ФНЧ в нули/полюса ФНЧ с другой частотой среза (в циклических частотах)</summary>
    /// <param name="Normed">Нормированные нули/полюса ФНЧ</param>
    /// <param name="wp">Новая циклическая частота среза</param>
    /// <returns>Нули/полюса нового ФНЧ</returns>
    public static IEnumerable<Complex> TransformToLowPassW(IEnumerable<Complex> Normed, double wp)
    {
        foreach (var p in Normed)
            yield return wp * p;
    }

    /// <summary>Преобразование нулей/полюсов нормированного ФНЧ в нули/полюса ФВЧ с другой частотой среза</summary>
    /// <param name="Normed">Нормированные нули/полюса ФНЧ</param>
    /// <param name="fp">Новая частота среза</param>
    /// <returns>Нули/полюса нового ФВЧ</returns>
    public static IEnumerable<Complex> TransformToHighPass(IEnumerable<Complex> Normed, double fp) =>
        TransformToHighPassW(Normed, pi2 * fp);

    /// <summary>Преобразование нулей/полюсов нормированного ФНЧ в нули/полюса ФВЧ с другой частотой среза (в циклических частотах)</summary>
    /// <param name="Normed">Нормированные нули/полюса ФНЧ</param>
    /// <param name="wp">Новая циклическая частота среза</param>
    /// <returns>Нули/полюса нового ФВЧ</returns>
    public static IEnumerable<Complex> TransformToHighPassW(IEnumerable<Complex> Normed, double wp)
    {
        foreach (var p in Normed)
            yield return wp / p;
    }
}