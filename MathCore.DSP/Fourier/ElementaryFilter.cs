using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathCore.DSP.Signals;

namespace MathCore.DSP.Fourier
{
    public class ElementaryFilter
    {
        private double _ReValue;
        private double _ImValue;
        private readonly int _N;
        private readonly int _m;
        private readonly double _dPhi;
        private int _n;

        public Complex Value => new Complex(_ReValue / _N, _ImValue / _N);

        public ElementaryFilter(int m, int N)
        {
            _m = m;
            _N = N;
            _dPhi = _m * Consts.pi2 / N;
        }

        public void Initialize()
        {
            _ReValue = 0;
            _ImValue = 0;
            _n = 0;
        }

        public void Process(double value)
        {
            var arg = _dPhi * _n++;
            _ReValue += Math.Cos(arg) * value;
            _ImValue += Math.Sin(arg) * value;
        }
    }
}
