﻿using System.Diagnostics;
using static System.Math;

namespace MathCore.DSP.Filters;

public class EllipticHighPass : EllipticFilter
{
    public static (double[] A, double[] B) GetPolynoms(Specification Spec)
    {
        var k_w = 1 / Spec.kw;
        var k_eps = 1 / Spec.kEps;

        var Kw = K(k_w);
        var Tw = T(k_w);
        var K_eps = K(k_eps);
        var T_eps = T(k_eps);

        var double_N = T_eps * Kw / (K_eps * Tw);
        var N = (int)Ceiling(double_N); // Порядок фильтра
        Debug.Assert(N > 0, $"N > 0 :: {N} > 0");

        var (zeros, poles) = GetNormedZerosPoles(N, Spec.EpsP, Spec.EpsS);

        // Масштабируем полюса на требуемую частоту пропускания
        var Ws = Spec.Ws;
        var high_pass_zeros_enum = TransformToHighPassW(zeros, Ws);
        if (N.IsOdd())
            high_pass_zeros_enum = high_pass_zeros_enum.AppendFirst(0);
        var high_pass_zeros = high_pass_zeros_enum;
        var high_pass_poles = TransformToHighPassW(poles, Ws);

        var zz = high_pass_zeros.ToArray();
        var pp = high_pass_poles.ToArray();

        // Переходим из p-плоскости в z-плоскость
        var dt = Spec.dt;
        var z_zeros = ToZArray(zz, dt);
        var z_poles = ToZArray(pp, dt);

        var k_zeros = z_zeros.Multiply(z => 1 + z).Re;
        var k_poles = z_poles.Multiply(z => 1 + z).Re;

        var g_norm = N.IsEven() 
            ? Spec.Gp * k_poles / k_zeros 
            : k_poles / k_zeros;


        var B = Polynom.Array.GetCoefficientsInverted(z_zeros).ToArray(v => v.Re * g_norm);
        var A = Polynom.Array.GetCoefficientsInverted(z_poles).ToRe();

        return (A, B);
    }

    public EllipticHighPass(
        double dt,
        double fs,
        double fp,
        double Gp = 0.891250938,
        double Gs = 0.01) 
        : this(GetSpecification(dt, fs, fp, Gp, Gs)) { }

    public EllipticHighPass(Specification Spec) : this(GetPolynoms(Spec), Spec) { }

    public EllipticHighPass((double[] A, double[] B) config, Specification Spec) : base(config.B, config.A, Spec) { }
}