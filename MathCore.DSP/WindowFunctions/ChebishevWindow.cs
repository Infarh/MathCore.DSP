using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathCore.DSP.WindowFunctions;

/// <summary>Параметрическое окно Дольфа-Чебышева</summary>
public static class ChebyshevWindow
{
    /// <summary>Значение отсчёта окна Чебышева</summary>
    /// <param name="n">Номер отсчёта окна (должен быть меньше N)</param>
    /// <param name="N">размер окна</param>
    /// <param name="gamma">Уровень боковых лепестков в дБ</param>
    /// <returns></returns>
    public static double Value(int n, int N, int gamma)
    {
        var q = 10.Pow(gamma / 20);
        
        var (m, d) = N % 2 == 0
            ? (N / 2 - 1, 0.5)  // N - чётное
            : ((N - 1) / 2, 0); // N - нечётное

#if NET8_0_OR_GREATER
        var b = Math.Cosh(Math.Acosh(q) / (N - 1));
#else
        var b = Math.Cosh(MathEx.Hyperbolic.Acosh(q) / (N - 1));
#endif

        static double C(int N1, double x) => Math.Cos(N1 * Math.Acos(x));

        var pin = Consts.pi / N;
        var p2in = Consts.pi2 / N;

        var sum = 0d;
        for(var i = 1; i <= m; i++)
        {
            var x1 = b * Math.Cos(pin * i);
            var x2 = Math.Cos(p2in * (n - m - d) * i);
            var c = C(N - 1, x1 * x2);
            sum += c;
        }

        var result = q + 2 * sum;

        return result;
    }
}
