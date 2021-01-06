using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
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
            /// <param name="Z">Нули/полюса прототипа (нормированного ФНЧ)</param>
            /// <param name="wpl">Нижняя частота среза ППФ</param>
            /// <param name="wph">Верхняя частота среза ППФ</param>
            /// <returns>Нули/полюса ППФ</returns>
            public static IEnumerable<Complex> ToBandPass(IEnumerable<Complex> Z, double wpl, double wph)
            {
                var dw = wph - wph;
                var wc2 = (wpl * wph);

                foreach (var z in Z)
                {
                    var dwz = 0.5 * dw * z;
                    var D = (dwz * dwz - wc2).Sqrt();

                    yield return -dwz + D;
                    yield return -dwz - D;
                }
            }
        }

        public readonly ref struct Specisication
        {
            public double dt { get; }

            public double fp { get; }
            public double fs { get; }

            public double wp { get; }
            public double ws { get; }

            public double Fp { get; }
            public double Fs { get; }

            public double Wp { get; }
            public double Ws { get; }

            public double Gp { get; }
            public double Gs { get; }

            public double Rp { get; }
            public double Rs { get; }

            public double EpsP { get; }
            public double EpsS { get; }

            public double kEps => EpsS / EpsP;
            public double kw => fs / fp;
            public double kW => Fs / Fp;

            public Specisication(double dt, double fp, double fs, double Gp, double Gs)
            {
                if (!(fp < fs)) throw new InvalidOperationException("Частота пропускания должна быть меньше частоты подавления");
                if (!(fp < 1 / (2 * dt))) throw new InvalidOperationException();

                this.dt = dt;
                this.fp = fp;
                this.fs = fs;
                this.Gp = Gp;
                this.Gs = Gs;

                Rp = -Gp.In_dB();
                Rs = -Gs.In_dB();

                EpsP = (10d.Pow(Rp / 10) - 1).Sqrt();
                EpsS = (10d.Pow(Rs / 10) - 1).Sqrt();

                Fp = ToAnalogFrequency(fp, dt);  // Частота пропускания аналогового прототипа
                Fs = ToAnalogFrequency(fs, dt);  // Частота подавления аналогового пропита

                wp = Consts.pi2 * fp;
                ws = Consts.pi2 * fp;

                Wp = Consts.pi2 * Fp;
                Ws = Consts.pi2 * Fp;
            }
        }

        public static Specisication GetSpecisication(double dt, double fp, double fs, double Gp, double Gs) => new(dt, fp, fs, Gp, Gs);

        /// <summary>Инициализация параметров цифрового фильтра на базе аналогового прототипа</summary>
        /// <param name="B">Коэффициенты полинома числителя</param>
        /// <param name="A">Коэффициенты полинома знаменателя</param>
        protected AnalogBasedFilter([NotNull] double[] B, [NotNull] double[] A) : base(B, A) { }
    }
}