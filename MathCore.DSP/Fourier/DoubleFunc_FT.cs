using System;

namespace MathCore.DSP.Fourier
{
    public static class DoubleFunc_FT
    {
        public static Func<double, Complex> GetFurierTransformation(this Func<double, double> s,
                                                                    double t1, double t2,
                                                                    bool IsInverse = false,
                                                                    double dt = 1e-4)
        {
            var lv_DeltaT = t2 - t1;
            var N = (int)Math.Abs(lv_DeltaT / dt);
            dt = lv_DeltaT / N;
            var w = Consts.pi2;
            if(IsInverse) w *= -1;

            return f =>
                       {
                           var pif = w * f;
                           if(IsInverse) pif *= -1;
                           var t = t1;
                           var val = s(t);
                           var arg = pif * t;
                           var p = val * Math.Cos(arg);
                           var q = val * Math.Sin(arg);

                           var reS = .0;
                           var imS = .0;
                           for(var i = 0; i < N; i++)
                           {
                               val = s(t);
                               arg = pif * t;
                               reS += p + (p = val * Math.Cos(arg));
                               imS += q + (q = val * Math.Sin(arg));
                               t += dt;
                           }
                           return new Complex(reS * .5 * dt, imS * .5 * dt);
                       };
        }

        public static Func<double, Complex> GetFurierTransformation(this Func<double, Complex> s,
                                                                    double t1, double t2,
                                                                    bool IsInverse = false,
                                                                    double dt = 1e-4)
        {
            var lv_DeltaT = t2 - t1;
            var N = (int)Math.Abs(lv_DeltaT / dt);
            dt = lv_DeltaT / N;
            var w = Consts.pi2;
            if(IsInverse) w *= -1;

            return f =>
                       {   //todo: проверить алгоритм!
                           var pif = w * f;
                           var t = t1;
                           var z = s(t).Rotate(pif * t);

                           var re = z.Re;
                           var im = z.Im;
                           for(var i = 0; i < N; i++)
                           {
                               z = s(t).Rotate(pif * t);
                               re += z.Re;
                               im += z.Im;
                               t += dt;
                           }
                           return (.5 * dt) * new Complex(re, im);
                       };
        }

        public static Func<int, Complex> GetFurierSpectrum(this Func<double, double> s,
                                                            double t1, double t2,
                                                            bool IsInverse = false,
                                                            double dt = 1e-4)
        {
            var lv_DeltaT = t2 - t1;
            var N = (int)Math.Abs(lv_DeltaT / dt);
            var ss = new double[N];
            dt = lv_DeltaT / N;
            for(var n = 0; n < N; n++)
                ss[n] = s(t1 + n * dt);
            return ss.GetFourierTransformation(IsInverse);


            //return m =>
            //           {
            //               var f = m / lv_DeltaT;
            //               var pif = Consts.pi2 * f;
            //               var t = t1;
            //               var val = s(t);
            //               var arg = pif * t;
            //               var p = val * Math.Cos(arg);
            //               var q = val * Math.Sin(arg);
            //               var reS = .0;
            //               var imS = .0;
            //               for(var i = 0; i < N; i++)
            //               {
            //                   val = s(t);
            //                   arg = pif * t;
            //                   reS += p + (p = val * Math.Cos(arg));
            //                   imS += q + (q = val * Math.Sin(arg));
            //                   t += dt;
            //               }
            //               return new Complex(reS * .5 * dt, imS * .5 * dt);
            //           };
        }

        public static Func<int, Complex> GetFurierSpectrum(this Func<double, Complex> s,
                                                            double t1, double t2,
                                                            bool IsInverse = false,
                                                            double dt = 1e-4)
        {
            var lv_DeltaT = t2 - t1;
            var N = (int)Math.Abs(lv_DeltaT / dt);

            var ss = new Complex[N];
            dt = lv_DeltaT / N;
            for(var n = 0; n < N; n++)
                ss[n] = s(t1 + n * dt);
            return ss.GetFourierTransformation(IsInverse);

            //return m =>
            //           {
            //               var f = m / lv_DeltaT;
            //               var pif = Consts.pi2 * f;
            //               var t = t1;
            //               var z = s(t).GetRotated(pif * t);
            //               var result = new Complex();
            //               for(var i = 0; i < N; i++)
            //               {
            //                   result.Inc(z);
            //                   result.Inc(z = s(t).GetRotated(pif * t));
            //                   t += dt;
            //               }
            //               return (.5 * dt) * result; //new Complex(reS * .5 * dt, imS * .5 * dt);
            //           };
        }
    }
}