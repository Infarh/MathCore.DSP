using System.Diagnostics.CodeAnalysis;
using System.Runtime.Serialization;

using static System.Math;
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters;

/// <summary>Фильтр Чебышева</summary>
[KnownType(typeof(ChebyshevLowPass))]
[KnownType(typeof(ChebyshevBandStop))]
public abstract class ChebyshevFilter : AnalogBasedFilter
{
    protected static double arcsh(double x) => Log(x + Sqrt(x * x + 1));
    protected static double arcch(double x) => Log(x + Sqrt(x * x - 1));

    /// <summary>Типы фильтров Чебышева</summary>
    public enum ChebyshevType : byte
    {
        /// <summary>Фильтр Чебышева первого рода - основной фильтр, пропускающий нижнюю полосу частот</summary>
        I,
        /// <summary>Фильтр Чебышева второго рода, подавляющий верхнюю область частот (выше fp)</summary>
        II,
        /// <summary>Фильтр Чебышева второго рода c коррекцией частотного диапазона, подавляющий верхнюю область частот (выше fs)</summary>
        IICorrected
    }

    protected static IEnumerable<Complex> GetNormedPolesI(int N, double EpsP, double W0 = 1)
    {
        var beta = arcsh(1 / EpsP) / N;
        var sh = Sinh(beta) * W0;
        var ch = Cosh(beta) * W0;
        if (N.IsOdd()) yield return -sh;    // Если порядок фильтра нечётный, то первым добавляем центральный полюс
        var r = N % 2;                   // Нечётность порядка фильтра
        for (var (i, dth) = (r, Consts.pi05 / N); i < N; i += 2)   // Расчёт полюсов
        {
            var th = dth * (i - r + 1);

            var sin = Sin(th);
            var cos = Cos(th);
            yield return new Complex(-sh * sin, +ch * cos);
            yield return new Complex(-sh * sin, -ch * cos);
        }
    }

    protected static (Complex[] Zeros, Complex[] Poles) GetNormedPolesII(int N, double EpsS, double W0 = 1)
    {
        var r = N % 2;                              // Нечётность порядка фильтра
        var L = N / 2;                              // Число пар нулей

        var beta = arcsh(EpsS) / N;
        var shb = Sinh(beta);
        var chb = Cosh(beta);

        var poles = new Complex[N];                 // Массив полюсов фильтра
        if (r != 0) poles[0] = -1 / shb;            // Если порядок фильтра нечётный, то первым добавляем центральный полюс
        for (var (i, dth) = (r, 0.5 * PI / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var th = dth * (i - r + 1);

            var sin = Sin(th);
            var cos = Cos(th);
            var norm = 1 / (sin * sin * shb * shb + cos * cos * chb * chb);
            //poles[i] = new Complex(-shb * sin * norm, chb * cos * norm);
            //poles[i + 1] = poles[i].ComplexConjugate;
            (poles[i], poles[i + 1]) = Complex.Conjugate(-shb * sin * norm, chb * cos * norm);
        }

        var zeros = new Complex[L * 2];
        for (var (n, dth) = (1, PI / N); n <= L; n++)
        {
            var th = dth * (n - 0.5);
            //zeros[2 * n - 2] = new Complex(0, 1 / Cos(th));
            //zeros[2 * n - 1] = zeros[2 * n - 2].ComplexConjugate;
            (zeros[2 * n - 2], zeros[2 * n - 1]) = Complex.Conjugate(0, 1 / Cos(th));
        }

        return (zeros, poles);
    }

    /// <summary>Тип фильтра</summary>
    public ChebyshevType FilterType { get; }

    /// <inheritdoc />
    /// <param name="Type">Тип фильтра I или II</param>
    protected ChebyshevFilter(double[] B, double[] A, Specification Spec, ChebyshevType Type) : base(B, A, Spec) => FilterType = Type;
}