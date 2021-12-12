using System.Collections;
using System.Collections.Generic;

namespace MathCore.DSP.Signals.IO
{
    public abstract class FileSignal : IEnumerable<double>
    {
        protected double _dt = double.NaN;
        protected double _t0;

        public string FileName { get; }

        public double dt
        {
            get
            {
                if (_dt > 0) return _dt;
                Initialize();
                return _dt;
            }
            set => _dt = value;
        }

        public double t0
        {
            get
            {
                if (_dt > 0) return _t0;
                Initialize();
                return _t0;
            }
            set => _t0 = value;
        }

        protected FileSignal(string FileName) => this.FileName = FileName;

        protected abstract void Initialize();

        public abstract IEnumerable<double> GetSamples();

        public abstract void SetSamples(IEnumerable<double> Samples);

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

        public IEnumerator<double> GetEnumerator() => GetSamples().GetEnumerator();
    }
}
