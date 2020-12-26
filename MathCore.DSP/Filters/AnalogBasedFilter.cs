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

        /// <summary>Инициализация параметров цифрового фильтра на базе аналогового прототипа</summary>
        /// <param name="B">Коэффициенты полинома числителя</param>
        /// <param name="A">Коэффициенты полинома знаменателя</param>
        protected AnalogBasedFilter([NotNull] double[] B, [NotNull] double[] A) : base(B, A) { }
    }
}