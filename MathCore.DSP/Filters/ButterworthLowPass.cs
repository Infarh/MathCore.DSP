using System;
using System.Linq;

namespace MathCore.DSP.Filters
{
    /// <summary>Низкочастотный фильтр Баттерворта</summary>
    public class ButterworthLowPass : ButterworthFilter
    {
        // <summary>Конфигурация низкочастотного фильтра Баттерворта, содержащая набор коэффициентов прямой и обратной связи</summary>
        // https://ru.dsplib.org/content/filter_butter_ap/filter_butter_ap.html

        /// <summary>Инициализация коэффициентов передаточной функции фильтра Баттерворта</summary>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="fp">Частота пропускания</param>
        /// <param name="fs">Частота заграждения</param>
        /// <param name="Gp">Затухание в полосе пропускания</param>
        /// <param name="Gs">Затухание в полосе заграждения</param>
        /// <returns>Кортеж с коэффициентами полинома числителя и знаменателя передаточной функции</returns>
        private static (double[] A, double[] B) Initialize(Specisication opt)
        {
            // Порядок фильтра
            var N = (int)Math.Ceiling(Math.Log(opt.kEps) / Math.Log(opt.kW));
            var r = N % 2; // Нечётность порядка фильтра

            // Радиус окружности размещения полюсов фильтра
            var alpha = opt.EpsP.Pow(-1d / N);

            // Угловой шаг между полюсами
            var dth = Math.PI / N;

            // Массив полюсов фильтра
            var poles = new Complex[N];
            // Если порядок фильтра нечётный, то первым добавляем центральный полюс
            if (r != 0) poles[0] = -alpha;
            // Расчёт полюсов
            for (var i = r; i < poles.Length; i += 2)
            {
                var w = dth * (i + 1 - r - 0.5);
                var sin = -alpha * Math.Sin(w);
                var cos = alpha * Math.Cos(w);
                poles[i] = new Complex(sin, cos);
                poles[i + 1] = new Complex(sin, -cos);
            }

            // Масштабируем полюса на требуемую частоту пропускания
            var Wp = opt.Wp;
            var translated_poles = poles.ToArray(p => p * Wp);
            // Переходим из p-плоскости в z-плоскость
            var dt = opt.dt;
            var z_poles = translated_poles.ToArray(p => ToZ(p, dt));
            // Вычисляем нормирующий множитель
            var kz = GetNormalizeCoefficient(translated_poles, dt);
            var WpN = Wp.Pow(N);
            var k = WpN * kz / opt.EpsP;
            var B = new double[N + 1];
            for (var i = 0; i < B.Length; i++)
                B[i] = k * SpecialFunctions.BinomialCoefficient(N, i);
            var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

            return (A, B);
        }

        /// <summary>Инициализация нового фильтра Баттерворта нижних частот</summary>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="fp">Частота пропускания</param>
        /// <param name="fs">Частота заграждения</param>
        /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
        /// <param name="Gs">Затухание в полосе заграждения (0.031622777 = -30 дБ)</param>
        public ButterworthLowPass(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.031622777)
            : this(Initialize(GetSpecisication(dt, fp, fs, Gp, Gs))) { }

        private ButterworthLowPass((double[] A, double[] B) config) : base(config.B, config.A) { }
    }
}
