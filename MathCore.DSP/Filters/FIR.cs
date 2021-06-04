using System;
using System.Collections.ObjectModel;

namespace MathCore.DSP.Filters
{
    /// <summary>Фильтр с конечной импульсной характеристикой</summary>
    public class FIR : DigitalFilter
    {
        /// <summary>Импульсная характеристика</summary>
        private readonly double[] _ImpulseResponse;

        /// <summary>Импульсная характеристика</summary>
        public ReadOnlyCollection<double> ImpulseResponse => Array.AsReadOnly(_ImpulseResponse);

        /// <summary>Инициализация нового цифрового фильтра с конечной импульсной характеристикой</summary>
        /// <param name="ImpulseResponse">Отсчёты импульсной характеристики фильтра</param>
        public FIR(double[] ImpulseResponse)
            : base((ImpulseResponse ?? throw new ArgumentNullException(nameof(ImpulseResponse))).Length) 
            => _ImpulseResponse = ImpulseResponse;

        public override double Process(double Sample, double[] state) => state.FilterSample(_ImpulseResponse, Sample);

        public override double Process(double Sample) => Process(Sample, State);

        public override Complex GetTransmissionCoefficient(double f) => _ImpulseResponse.GetTransmissionCoefficient(f);
    }
}