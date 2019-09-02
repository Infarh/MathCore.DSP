using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Linq;
using MathCore.Annotations;

namespace MathCore.DSP
{
    public static class DoubleArrayDSPExtensions
    {
        public static Complex GetTransmissionCoefficient([NotNull] this double[] ImpulseResponse, double f, double dt)
            => ImpulseResponse.GetTransmissionCoefficient(f * dt);
        public static Complex GetTransmissionCoefficient([NotNull] this double[] ImpulseResponse, double f)
        {
            var e = Complex.Exp(-2 * Math.PI * f);
            Complex result = ImpulseResponse[ImpulseResponse.Length - 1];
            for (var i = ImpulseResponse.Length - 2; i >= 0; i--)
                result = result * e + ImpulseResponse[i];
            return result;
        }

        public static double FilterSample([NotNull] this double[] State, [NotNull] double[] ImpulseResponse, double Sample)
        {
            Contract.Ensures(State != null);
            Contract.Ensures(ImpulseResponse != null);
            Contract.Ensures(State.Length == ImpulseResponse.Length);

            var result = 0d;

            for (var i = State.Length - 1; i >= 1; i--)
            {
                State[i] = State[i - 1];
                result += State[i] * ImpulseResponse[i];
            }

            State[0] = Sample;

            return result + Sample * ImpulseResponse[0];
        }

        [NotNull]
        public static IEnumerable<double> FilterFIR([NotNull] this IEnumerable<double> samples, [NotNull] double[] ImpulseResponse, [NotNull] double[] State)
        {
            if (samples is null) throw new ArgumentNullException(nameof(samples));
            if (ImpulseResponse is null) throw new ArgumentNullException(nameof(ImpulseResponse));
            if (State is null) throw new ArgumentNullException(nameof(State));
            if (ImpulseResponse.Length != State.Length) throw new InvalidOperationException("Размер массива импульсной характеристики не соответствует размеру массива состояния фильтра");
            return samples.Select(sample => State.FilterSample(ImpulseResponse, sample));
        }

        [NotNull]
        public static IEnumerable<double> FilterFIR([NotNull] this IEnumerable<double> samples, [NotNull] double[] ImpulseResponse)
            => samples.FilterFIR(ImpulseResponse, new double[ImpulseResponse.Length]);


        public static Complex GetTransmissionCoefficient([NotNull] double[] A, [NotNull] double[] B, double f, double dt)
            => GetTransmissionCoefficient(A, B, f * dt);

        public static Complex GetTransmissionCoefficient([NotNull] double[]A, [NotNull] double[]B, double f)
        {
            var e = Complex.Exp(-2 * Math.PI * f);

            Complex Sum(double[] V, Complex exp)
            {
                Complex sum = V[V.Length - 1];
                for (var i = V.Length - 2; i >= 0; i--) sum = sum * exp + V[i];
                return sum;
            }

            return Sum(B, e) / Sum(A, e);
        }

        public static double FilterSample([NotNull] this double[] State, [NotNull] double[] A, [NotNull] double[] B, double Sample)
        {
            Contract.Requires(State != null);
            Contract.Requires(A != null);
            Contract.Requires(B != null);
            Contract.Requires(A.Length >= B.Length);
            Contract.Requires(A.Length > 2);
            Contract.Requires(B.Length > 1);
            Contract.Requires(A[0] != 0);

            var a0 = 1 / A[0];

            var result = 0d;
            var input = Sample;
            var b_length = B.Length;
            if (A.Length == b_length)
                for (var i = State.Length - 1; i >= 1; i--)
                {
                    var v = State[i - 1];
                    State[i] = v;
                    result += v * B[i] * a0;
                    input -= v * A[i] * a0;
                }
            else
            {
                for (var i = State.Length - 1; i >= b_length; i--)
                {
                    var v = State[i - 1];
                    State[i] = v;
                    input -= v * A[i] * a0;
                }
                for (var i = b_length - 1; i >= 1; i--)
                {
                    var v = State[i - 1];
                    State[i] = v;
                    result += v * B[i] * a0;
                    input -= v * A[i] * a0;
                }
            }

            State[0] = input;
            var r = result + input * B[0] * a0;
            return r;
        }

        [NotNull]
        public static IEnumerable<double> FilterIIR([NotNull] this IEnumerable<double> samples, [NotNull] double[] A, [NotNull] double[] B, [NotNull] double[] State)
        {
            if (samples is null) throw new ArgumentNullException(nameof(samples));
            if (A is null) throw new ArgumentNullException(nameof(A));
            if (B is null) throw new ArgumentNullException(nameof(B));
            if (A.Length < B.Length) throw new InvalidOperationException("Размеры массивов числителя и знаменателя передаточной функции не равны");
            if (State is null) throw new ArgumentNullException(nameof(State));
            return samples.Select(sample => FilterSample(State, A, B, sample));
        }

        [NotNull] public static IEnumerable<double> FilterIIR(this IEnumerable<double> samples, [NotNull] double[] A, [NotNull] double[] B)
            => samples.FilterIIR(A, B, new double[A.Length]);
    }
}
