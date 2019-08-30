using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MathCore.DSP.Fourier;
using MathCore.DSP.Signals.Operations;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Signals
{
    public abstract class DigitalSignal : IEnumerable<double>
    {
        protected readonly double _dt;
        public double dt => _dt;
        public double fd => 1 / _dt;

        public double TotalTime => SamplesCount * _dt;

        public abstract int SamplesCount { get; }

        public double Min => this.Min();
        public double Max => this.Max();

        public double PeakToPeakAmplitude => this.GetMinMax().Length;

        public double Power => this.Average(s => s * s);

        public abstract double this[int n] { get; set; }

        protected DigitalSignal(double dt)
        {
            if (dt <= 0) throw new ArgumentOutOfRangeException(nameof(dt), "Период дискретизации должен быть больше 0");
            _dt = dt;
        }

        public abstract IEnumerator<double> GetEnumerator();
        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

        #region Операторы

        public static SumOfSignalsResultSignal operator +(DigitalSignal s1, DigitalSignal s2) => new SumOfSignalsResultSignal(s1, s2);
        public static SubstractionOfSignalsResultSignal operator -(DigitalSignal s1, DigitalSignal s2) => new SubstractionOfSignalsResultSignal(s1, s2);
        public static MultiplyOfSignalsResultSignal operator *(DigitalSignal s1, DigitalSignal s2) => new MultiplyOfSignalsResultSignal(s1, s2);
        public static DivisionOfSignalsResultSignal operator /(DigitalSignal s1, DigitalSignal s2) => new DivisionOfSignalsResultSignal(s1, s2);

        public static SumOfSignalWithScalarResultSignal operator +(DigitalSignal s, double x) => new SumOfSignalWithScalarResultSignal(s, x);
        public static SumOfSignalWithScalarResultSignal operator +(double x, DigitalSignal s) => new SumOfSignalWithScalarResultSignal(s, x);
        public static SubstractionOfSignalWithScalarResultSignal operator -(DigitalSignal s, double x) => new SubstractionOfSignalWithScalarResultSignal(s, x);
        public static SubstractionOfScalarWithSignalResultSignal operator -(double x, DigitalSignal s) => new SubstractionOfScalarWithSignalResultSignal(s, x);
        public static MultiplyOfSignalWithScalarResultSignal operator *(DigitalSignal s, double x) => new MultiplyOfSignalWithScalarResultSignal(s, x);
        public static MultiplyOfSignalWithScalarResultSignal operator *(double x, DigitalSignal s) => new MultiplyOfSignalWithScalarResultSignal(s, x);
        public static DivisionOfSignalWithScalarResultSignal operator /(DigitalSignal s, double x) => new DivisionOfSignalWithScalarResultSignal(s, x);
        public static DivisionOfScalarWithSignalResultSignal operator /(double x, DigitalSignal s) => new DivisionOfScalarWithSignalResultSignal(s, x);

        #endregion

        //public virtual Complex[] GetSpectrumSamples(double dt) => FT.FourierTransform()
    }

}