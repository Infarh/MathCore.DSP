using System.Runtime.Serialization;

using static MathCore.Consts;
// ReSharper disable InconsistentNaming
// ReSharper disable ParameterHidesMember

namespace MathCore.DSP.Filters;

/// <summary>Цифровой фильтр на основе аналогового прототипа</summary>
//[KnownType(typeof(BesselFilter))]       // Наихудшее аппроксимация прямоугольной АЧХ (близка по форме к гаусовой кривой), но наилучшая форма переходной хар-ки
[KnownType(typeof(ButterworthFilter))]  // Промежуточные качества по прямоугольности АЧХ и переходной хар-ке
[KnownType(typeof(ChebyshevFilter))]    // АЧХ наиболее приближена к прямоугольной - наибольшие пульсации в переходной хор-ке
[KnownType(typeof(EllipticFilter))]     // Максимальная крутизна АЧХ
public abstract class AnalogBasedFilter : IIR
{
    //https://ru.dsplib.org/content/filter_low2low/filter_low2low.html
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
        public static IEnumerable<Complex> ToHigh(IEnumerable<Complex> Z, double wp) => Z.Select(z => wp * z);

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

            var wpl = pi2 * fpl;
            var wph = pi2 * fph;

            var dw = (wph - wph) / 2;
            var wc2 = wpl * wph;

            // На каждый исходный полюс формируется пара новых полюсов
            // Число нулей равно удвоенному числу исходных нулей
            //      + ноль в 0 кратности, равной разности числа полюсов и нулей
            var poles = new Complex[Poles.Length * 2];
            var zeros = new Complex[Zeros.Length * 2 + (Poles.Length - Zeros.Length)];

            static void Set(in Complex p, in Complex D, out Complex p1, out Complex p2)
            {
                p1 = p + D;
                p2 = p - D;
            }

            for (var i = 0; i < Poles.Length; i++)
            {
                var pdw = dw * Poles[i];
                Set(pdw, Complex.Sqrt(pdw.Power - wc2),
                    out poles[2 * i],
                    out poles[2 * i + 1]);
            }

            for (var i = 0; i < Zeros.Length; i++)
            {
                var pdw = dw * Zeros[i];
                Set(pdw, Complex.Sqrt(pdw.Power - wc2),
                    out zeros[2 * i],
                    out zeros[2 * i + 1]);
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
        public double kW => Fs / Fp;

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
            if (!(fp < fs)) throw new InvalidOperationException("Частота пропускания должна быть меньше частоты подавления");
            if (!(fp < 1 / (2 * dt))) throw new InvalidOperationException("Частота пропускания должен быть меньше половины полосы дискретизации");

            this.dt = dt;
            this.fp = fp;
            this.fs = fs;
            this.Gp = Gp;
            this.Gs = Gs;

            Rp = -Gp.In_dB();
            Rs = -Gs.In_dB();

            EpsP = (10d.Pow(Rp / 10) - 1).Sqrt();
            EpsS = (10d.Pow(Rs / 10) - 1).Sqrt();

            //var tEpsP = (1 / Gp.Pow2() - 1).Sqrt();
            //var tEpsS = (1 / Gs.Pow2() - 1).Sqrt();

            //var tEpsP = (1 - Gp*Gp).Sqrt() / (Gp*Gp);
            //var tEpsS = (1 - Gs*Gs).Sqrt() / (Gs*Gs);

            Fp = ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
            Fs = ToAnalogFrequency(fs, dt);  // Частота подавления аналогового пропита

            wp = pi2 * fp;
            ws = pi2 * fp;

            Wp = pi2 * Fp;
            Ws = pi2 * Fp;
        }

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

    public override Complex GetTransmissionCoefficient(double f) => base.GetTransmissionCoefficient(f / fd);

    public override Complex GetTransmissionCoefficient(double f, double dt) => base.GetTransmissionCoefficient(f / fd);

    public static IEnumerable<Complex> TransformToBandPassPoles(IEnumerable<Complex> NormedPoles, double fmin, double fmax)
    {
        var w_min = Consts.pi2 * fmin;
        var w_max = Consts.pi2 * fmax;
        var dw = (w_max - w_min) / 2;
        var w2 = w_min * w_max;

        foreach (var p in NormedPoles)
        {
            var pdw = p * dw;
            var sqrt = Complex.Sqrt(pdw.Pow2() - w2);
            yield return pdw + sqrt;
            yield return pdw - sqrt;
        }
    }
}