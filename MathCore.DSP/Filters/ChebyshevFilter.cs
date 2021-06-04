using System.Diagnostics.CodeAnalysis;
using System.Runtime.Serialization;

using static System.Math;

namespace MathCore.DSP.Filters
{
    /// <summary>Фильтр Чебышева</summary>
    [KnownType(typeof(ChebyshevLowPass))]
    public abstract class ChebyshevFilter : AnalogBasedFilter
    {
        protected static double arcsh(double x) => Log(x + Sqrt(x * x + 1));
        protected static double arcch(double x) => Log(x + Sqrt(x * x - 1));

        /// <summary>Типы фильтров Чебышева</summary>
        [SuppressMessage("ReSharper", "InconsistentNaming")]
        public enum ChebyshevType : byte
        {
            /// <summary>Фильтр Чебышева первого рода - основной фильтр, пропускающий нижнюю полосу частот</summary>
            I,
            /// <summary>Фильтр Чебышева второго рода, подавляющий верхнюю область частот (выше fp)</summary>
            II,
            /// <summary>Фильтр Чебышева второго рода c коррекцией частотного диапазона, подавляющий верхнюю область частот (выше fs)</summary>
            IICorrected
        }

        protected static Complex[] GetNormedPolesI(int N, double EpsP)
        {
            var r = N % 2;                              // Нечётность порядка фильтра
            var dth = PI / N;                      // Угловой шаг между полюсами
            var beta = arcsh(1 / EpsP) / N;
            var sh = Sinh(beta);
            var ch = Cosh(beta);
            var poles = new Complex[N];                 // Массив полюсов фильтра
            if (r != 0) poles[0] = -sh;                 // Если порядок фильтра нечётный, то первым добавляем центральный полюс
            for (var i = r; i < poles.Length; i += 2)   // Расчёт полюсов
            {
                var n = (i - r) / 2 + 1;
                var th = dth * (n - 0.5);

                var sin = Sin(th);
                var cos = Cos(th);
                poles[i] = new Complex(-sh * sin, ch * cos);
                poles[i + 1] = poles[i].ComplexConjugate;
            }

            return poles;
        }

        protected static (Complex[] Zeros, Complex[] Poles) GetNormedPolesII(int N, double EpsS)
        {
            var r = N % 2;                              // Нечётность порядка фильтра
            var L = (N - r) / 2;                        // Число пар нулей
            var dth = PI / N;                      // Угловой шаг между полюсами
            var beta = arcsh(EpsS) / N;
            var shb = Sinh(beta);
            var chb = Cosh(beta);

            var poles = new Complex[N];                 // Массив полюсов фильтра
            if (r != 0) poles[0] = -1 / shb;            // Если порядок фильтра нечётный, то первым добавляем центральный полюс
            for (var i = r; i < poles.Length; i += 2)   // Расчёт полюсов
            {
                var n = (i - r) / 2 + 1;
                var th = dth * (n - 0.5);

                var sin = Sin(th);
                var cos = Cos(th);
                var norm = 1 / (sin * sin * shb * shb + cos * cos * chb * chb);
                poles[i] = new Complex(-shb * sin * norm, chb * cos * norm);
                poles[i + 1] = poles[i].ComplexConjugate;
            }

            var zeros = new Complex[L * 2];
            for (var n = 1; n <= L; n++)
            {
                var th = dth * (n - 0.5);
                zeros[2 * n - 2] = new Complex(0, 1 / Cos(th));
                zeros[2 * n - 1] = zeros[2 * n - 2].ComplexConjugate;
            }

            return (zeros, poles);
        }

        /// <inheritdoc />
        protected ChebyshevFilter(double[] B, double[] A) : base(B, A) { }
    }
}