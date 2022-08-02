using static System.Math;

using static MathCore.Polynom.Array;
using static MathCore.SpecialFunctions;

namespace MathCore.DSP.Filters;

/// <summary>Низкочастотный фильтр Баттерворта</summary>
public class ButterworthLowPass : ButterworthFilter
{
    public static int GetOrder(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.01)
    {
        var kEps2 = (1 / (Gs * Gs) - 1) / (1 / (Gp * Gp) - 1);
        var pi_dt = PI * dt;
        var kW = Tan(fs * pi_dt) / Tan(fp * pi_dt);

        var N = (int)Ceiling(0.5 * Log(kEps2) / Log(kW));
        return N;
    }

    public static double GetFrequencyStop(double dt, double fp, int Order, double Gp = 0.891250938, double Gs = 0.01)
    {
        var kEps2 = (1 / (Gs * Gs) - 1) / (1 / (Gp * Gp) - 1);
        var pi_dt = PI * dt;

        var fs = Atan(Pow(kEps2, 0.5 / Order) * Tan(pi_dt * fp)) / pi_dt;
        return fs;
    }

    public static double GetFrequencyPass(double dt, int Order, double fs, double Gp = 0.891250938, double Gs = 0.01)
    {
        var kEps2 = (1 / (Gs * Gs) - 1) / (1 / (Gp * Gp) - 1);
        var pi_dt = PI * dt;

        var fp = Atan(Pow(kEps2, -0.5 / Order) * Tan(pi_dt * fs)) / pi_dt;
        return fp;
    }

    public static IEnumerable<Complex> GetPoles(int Order, double Gp = 0.891250938) => GetNormPolesGp(Order, Gp);

    // https://ru.dsplib.org/content/filter_butter_ap/filter_butter_ap.html

    /// <summary>Инициализация коэффициентов передаточной функции фильтра Баттерворта</summary>
    /// <returns>Кортеж с коэффициентами полинома числителя и знаменателя передаточной функции</returns>
    public static (double[] A, double[] B) GetPolynoms(Specification Spec)
    {
        // Порядок фильтра
        var N = (int)Ceiling(Log(Spec.kEps) / Log(Spec.kW));
        var poles = GetNormPoles(N, Spec.EpsP);

        // Масштабируем полюса на требуемую частоту пропускания
        var Wp = Spec.Wp;
        var translated_poles = TransformToLowPassW(poles, Wp);

        // Переходим из p-плоскости в z-плоскость
        var dt = Spec.dt;
        var z_poles = ToZArray(translated_poles, dt);

        // Вычисляем нормирующий множитель
        var g_norm = z_poles.Multiply(z => (1 - z) / 2).Re;

        var B = new double[N + 1];
        for (var i = 0; i < B.Length; i++)
            B[i] = g_norm * BinomialCoefficient(N, i);

        var A = GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    /// <summary>Инициализация нового фильтра Баттерворта нижних частот</summary>
    /// <param name="dt">Период дискретизации</param>
    /// <param name="fp">Частота пропускания</param>
    /// <param name="fs">Частота заграждения</param>
    /// <param name="Gp">Затухание в полосе пропускания (0.891250938 = -1 дБ)</param>
    /// <param name="Gs">Затухание в полосе заграждения (0.01 = -40 дБ)</param>
    public ButterworthLowPass(double dt, double fp, double fs, double Gp = 0.891250938, double Gs = 0.01)
        : this(GetSpecification(dt, fp, fs, Gp, Gs)) { }

    /// <summary>Инициализация нового фильтра Баттерворта нижних частот</summary>
    /// <param name="Spec">Спецификация фильтра</param>
    public ButterworthLowPass(Specification Spec) : this(GetPolynoms(Spec), Spec) { }

    private ButterworthLowPass((double[] A, double[] B) config, Specification Spec) : base(config.B, config.A, Spec) { }
}