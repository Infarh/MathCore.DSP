using System;
using System.Collections.ObjectModel;
using MathCore.Annotations;

namespace MathCore.DSP.Filters
{
    /// <summary>Фильтр с бесконечной импульсной характеристикой</summary>
    public class IIR : DigitalFilter
    {
        /// <summary>Массив кооэффициентов полинома числителя</summary>
        [NotNull] private readonly double[] _B;
        /// <summary>Массив коэффициентов полинома знаменателя</summary>
        [NotNull] private readonly double[] _A;

        /// <summary>Массив кооэффициентов полинома числителя</summary>
        public ReadOnlyCollection<double> B => Array.AsReadOnly(_B);

        /// <summary>Массив коэффициентов полинома знаменателя</summary>
        public ReadOnlyCollection<double> A => Array.AsReadOnly(_A);

        /// <summary>Инициализация нового цифрового фильтра с бесконечной импульсной характеристикой</summary>
        /// <param name="B">Массив кооэффициентов полинома числителя</param>
        /// <param name="A">Массив коэффициентов полинома знаменателя</param>
        public IIR([NotNull] double[] B, [NotNull] double[] A)
            : base(
                Math.Max(
                    (A ?? throw new ArgumentNullException(nameof(B))).Length,
                    (B ?? throw new ArgumentNullException(nameof(A))).Length))
        {
            if (B.Length == 0) throw new ArgumentException("Размер массива коэффициентов числителя должен быть больше 0", nameof(B));
            if (A.Length < 2) throw new ArgumentException("Размер массива коэффициентов знаменателя должен быть больше 1", nameof(A));
            if (B.Length > A.Length) throw new ArgumentException("Размер массива коэффициентов полинома числителя должен быть меньше, либо равен размеру массива коэффициентов полинома знаменателя");
            _B = B;
            _A = A;
        }

        public override double Process(double Sample, double[] state) => state.FilterSample(_A, _B, Sample);

        public override double Process(double Sample) => Process(Sample, State);

        public override Complex GetTransmissionCoefficient(double f) => DoubleArrayDSPExtensions.GetTransmissionCoefficient(_A, _B, f);

        /// <summary>Последовательное соединение фильтра с другим <see cref="IIR"/></summary>
        /// <param name="filter">Соединяемый фильтр</param>
        /// <returns>Фильтр, представляющий собой результат последовательного соединения двух фильтров</returns>
        public IIR ConnectionSerialTo([NotNull] IIR filter)
        {
            if (filter is null) throw new ArgumentNullException(nameof(filter));

            var b = Polynom.Array.Multiply(_B, filter._B);
            var a = Polynom.Array.Multiply(_A, filter._A);
            return new IIR(b, a);
        }

        /// <summary>Параллельное соединение фильтра с другим <see cref="IIR"/></summary>
        /// <param name="filter">Соединяемый фильтр</param>
        /// <returns>Фильтр, представляющий собой результат параллельного соединения двух фильтров</returns>
        public IIR ConnectionParallelTo([NotNull] IIR filter)
        {
            if (filter is null) throw new ArgumentNullException(nameof(filter));

            var b1 = _B;
            var b2 = filter._B;

            var a1 = _A;
            var a2 = filter._A;

            var b = Polynom.Array.Sum(Polynom.Array.Multiply(b1, a2), Polynom.Array.Multiply(b2, a1));
            var a = Polynom.Array.Multiply(a1, a2);
            return new IIR(Polynom.Array.Multiply(b, 0.5), a);
        }
    }
}