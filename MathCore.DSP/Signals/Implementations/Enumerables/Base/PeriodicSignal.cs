// ReSharper disable InconsistentNaming
namespace MathCore.DSP.Signals.Implementations.Enumerables.Base;

/// <summary>Гармонический сигнал</summary>
public abstract class PeriodicSignal : EnumerableSignalImplementation
{
    /// <summary>Гармоническая основа</summary>
    protected abstract class PeriodicSignalInfo : SignalInfo
    {
        /// <summary>Частота</summary>
        protected double _f0;
        /// <summary>Фаза</summary>
        protected double _phi0;

        /// <summary>Частота</summary>
        public double f0 { get => _f0; set => _f0 = value; }

        /// <summary>Циклическая частота</summary>
        public double w0 { get => Consts.pi2 * _f0; set => _f0 = value / Consts.pi2; }

        /// <summary>Период</summary>
        public double T0 { get => 1 / _f0; set => _f0 = 1 / value; }

        /// <summary>Фаза</summary>
        public double phi0 { get => _phi0; set => _phi0 = value; }

        /// <summary>Временное смещение</summary>
        public double DeltaT { get => _phi0 / w0; set => _phi0 = value * w0; }

        /// <summary>Инициализация новой гармонической основы сигнала</summary>
        /// <param name="A">Амплитуда сигнала</param>
        /// <param name="f0">Частота</param>
        /// <param name="phi0">Начальная фаза</param>
        /// <param name="dt">Период дискретизации</param>
        /// <param name="SamplesCount">Число отсчётов сигнала</param>
        protected PeriodicSignalInfo(double A, double f0, double phi0, double dt, int SamplesCount) : base(A, dt, SamplesCount)
        {
            _f0 = f0;
            _phi0 = phi0;
        }
    }

    private readonly PeriodicSignalInfo _SignalInfo;

    public double A { get => _SignalInfo.A; set => _SignalInfo.A = value; }

    public double f0 { get => _SignalInfo.f0; set => _SignalInfo.f0 = value; }

    public double phi0 { get => _SignalInfo.phi0; set => _SignalInfo.phi0 = value; }
            
    protected PeriodicSignal(double dt, PeriodicSignalInfo SignalInfo, double t0) : base(dt, SignalInfo, t0) => _SignalInfo = SignalInfo;
}