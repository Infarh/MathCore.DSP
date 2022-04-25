using System.Diagnostics.CodeAnalysis;
using System.Runtime.Serialization;

using MathCore.DSP.Infrastructure;

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
        var beta = arcsh(EpsS) / N;
        var sh = Sinh(beta);
        var ch = Cosh(beta);

        var poles = new Complex[N];                 // Массив полюсов фильтра
        if (N.IsOdd())
            poles[0] = -W0 / sh;                    // Если порядок фильтра нечётный, то первым добавляем центральный полюс
        var r = N % 2;                              // Нечётность порядка фильтра
        for (var (i, dth) = (r, 0.5 * PI / N); i < poles.Length; i += 2)   // Расчёт полюсов
        {
            var (cos, sin) = Complex.Exp(dth * (i - r + 1));
            var z = new Complex(-sh * sin, ch * cos);
            var norm = W0 / z.Power;
            (poles[i], poles[i + 1]) = Complex.Conjugate(z.Re * norm, z.Im * norm);
        }

        var zeros = new Complex[N - r];
        for (var (n, dth, L) = (1, PI / N, N / 2); n <= L; n++)
            (zeros[2 * n - 2], zeros[2 * n - 1]) = Complex.Conjugate(0, W0 / Cos(dth * (n - 0.5)));

        return (zeros, poles);
    }

    /// <summary>Тип фильтра</summary>
    public ChebyshevType FilterType { get; }

    /// <inheritdoc />
    /// <param name="Type">Тип фильтра I или II</param>
    protected ChebyshevFilter(double[] B, double[] A, Specification Spec, ChebyshevType Type) : base(B, A, Spec) => FilterType = Type;
}