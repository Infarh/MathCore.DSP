
using System;

namespace MathCore.DSP.Fourier
{
    internal struct exp
    {
        public readonly double sin;

        public readonly double cos;

        public exp(double sin, double cos) { this.sin = sin; this.cos = cos; }

        public static exp[] GetCoefficients(int N, bool IsInverse = false)
        {
            var w = new exp[N];
            var darg = Consts.pi2 / N;
            if(IsInverse) darg *= -1;
            var arg = 0.0;
            for(var i = 0; i < N; i++)
            {
                w[i] = new exp(Math.Sin(arg), Math.Cos(arg));
                arg += darg;
            }
            return w;
        }
    }
}