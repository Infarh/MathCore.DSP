using System.Diagnostics.CodeAnalysis;

using static System.Math;

#pragma warning disable 1591

namespace MathCore.DSP.Fourier;

[Copyright("alglib", url = "alglib.sources.ru")]
#pragma warning disable CS8981 // Стили именования
#pragma warning disable IDE1006 // Стили именования
public static class fft
#pragma warning restore IDE1006 // Стили именования
#pragma warning restore CS8981 // Стили именования
{
    /// <summary>
    /// Быстрое комплексное одномерное преобразование Фурье
    /// </summary>
    /// <param name="x">Вектор комплексных значений, преобразуемый в спектр</param>
    /// <remarks>
    /// Размер массива N может быть произвольным числом (составным или простым). 
    /// Составной N будут обработаны с использованием вариативного алгоритма Кули-Тьюки с кешированием.
    /// Массивы, размер которых соответствует малому простому числу преобразуются с помощью жестко сода
    /// (по аналогии с кодом FFTW, но без оптимизации низкого уровня), большое простое число элементов 
    /// обрабатывается с помощью алгоритма Блустейна.
    /// 
    /// Быстрые преобразования для гладких N (только простые множители 2, 3, 5), самый быстрый для степеней 2. 
    /// При N имеющих простые множители большие, чем эти, но порядка меньше, чем N, вычисления 
    /// будут примерно в 4 раза медленнее, чем для близких высоко составных N. 
    /// Когда N является простым, скорость будет в 6 раз меньше.
    /// 
    /// Алгоритмическая сложность O(N*logN) для любых N
    /// </remarks>
    [Copyright("29.05.2009 by Bochkanov Sergey", url = "alglib.sources.ru")]
    public static Complex[] FFT(Complex[] x)
    {
        var N = x.Length;
        if (N == 1) return [x[0]];

        var buf = new double[2 * N];
        for (var i = 0; i < N; i++)
            (buf[2 * i], buf[2 * i + 1]) = x[i];

        //
        // Generate plan and execute it.
        //
        // Plan is a combination of a successive factorizations of N and
        // precomputed data. It is much like a FFTW plan, but is not stored
        // between subroutine calls and is much simpler.
        //
        var plan = new FT_Base.FT_Plan();
        FT_Base.FTBaseGenerateComplexFFTPlan(N, plan);
        FT_Base.FTBaseExecutePlan(ref buf, 0, plan);

        var result = new Complex[N];
        var last = 0d;
        for (int i = 0, N2 = buf.Length, i0 = 0; i < N2; i++)
            if (i % 2 == 0)
                last = buf[i] / N;
            else
                result[i0++] = new(last, buf[i] / N);
        return result;
    }


    /// <summary>
    /// Быстрое обратное комплексное одномерное преобразование Фурье
    /// </summary>
    /// <param name="y">Массив значений спектра</param>
    [Copyright("29.05.2009 by Bochkanov Sergey", url = "alglib.sources.ru")]
    public static Complex[] FFT_Complex_Inverse(Complex[] y)
    {
        var N = y.Length;

        var result = new Complex[N];
        for (var i = 0; i < N; i++)
            result[i] = y[i].ComplexConjugate;

        result = FFT(result);//, n);

        for (var i = 0; i < N; i++)
            result[i] = result[i].ComplexConjugate * N; // / N;
        return result;
    }


    /*************************************************************************
    1-dimensional real FFT.

    Algorithm has O(N*logN) complexity for any N (composite or prime).

    INPUT PARAMETERS
        A   -   array[0..N-1] - real function to be transformed
        N   -   problem size

    OUTPUT PARAMETERS
        F   -   DFT of a input array, array[0..N-1]
                F[j] = SUM(A[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)

    NOTE
        F[] satisfies symmetry property F[k] = conj(F[N-k]),  so just one half
    of  array  is  usually needed. But for convinience subroutine returns full
    Complex array (with frequencies above N/2), so its result may be  used  by
    other FFT-related subroutines.


      -- ALGLIB --
         Copyright 01.06.2009 by Bochkanov Sergey
    *************************************************************************/
    /// <summary>Быстрое одномерное вещественное преобразование Фурье</summary>
    /// <param name="x">Массив входных значений</param>
    /// <value>Массив комплексных значений спектра</value>
    [SuppressMessage("ReSharper", "TooWideLocalVariableScope")]
    public static Complex[] FFT(double[] x)
    {
        var N = x.Length;
        switch (x)
        {
            case [var x0]: return [new(x0)];
            case [var x0, var x1]: return [new(x0 + x1), new(x0 - x1)];
        }

        var result = new Complex[N];
        if (N % 2 == 0)
        {
            var plan = new FT_Base.FT_Plan();

            var n05 = N / 2;
            FT_Base.FTBaseGenerateComplexFFTPlan(n05, plan);
            var z = (double[])x.Clone();
            FT_Base.FTBaseExecutePlan(ref z, 0, plan);

            var pi_n = -Consts.pi2 / N;
            var N05 = .5 / N;
            int i_dx;
            double h_n_re, h_n_im;
            double h_mn_c_re, h_mn_c_im;
            double nsin, cos;
            for (var i = 0; i <= n05; i++)
            {
                nsin = -Sin(pi_n * i);
                cos = Cos(pi_n * i);

                i_dx = 2 * (i % n05);
                h_n_re = z[i_dx];
                h_n_im = z[i_dx + 1];

                i_dx = 2 * ((n05 - i) % n05);
                h_mn_c_re = z[i_dx];
                h_mn_c_im = -z[i_dx + 1];

                result[i] = new
                (
                    Re: (h_n_re + h_mn_c_re - nsin * (h_n_re - h_mn_c_re) + cos * (h_n_im - h_mn_c_im)) * N05,
                    Im: (h_n_im + h_mn_c_im - nsin * (h_n_im - h_mn_c_im) - cos * (h_n_re - h_mn_c_re)) * N05
                );
            }
            for (var i = n05 + 1; i < N; i++)
                result[i] = result[N - i].ComplexConjugate;
            return result;
        }

        for (var i = 0; i < N; i++)
            result[i] = new(x[i]);
        return FFT(result);
    }


    /*************************************************************************
    1-dimensional real inverse FFT.

    Algorithm has O(N*logN) complexity for any N (composite or prime).

    INPUT PARAMETERS
        F   -   array[0..floor(N/2)] - frequencies from forward real FFT
        N   -   problem size

    OUTPUT PARAMETERS
        A   -   inverse DFT of a input array, array[0..N-1]

    NOTE
        F[] should satisfy symmetry property F[k] = conj(F[N-k]), so just  one
    half of frequencies array is needed - elements from 0 to floor(N/2).  F[0]
    is ALWAYS real. If N is even F[floor(N/2)] is real too. If N is odd,  then
    F[floor(N/2)] has no special properties.

    Relying on properties noted above, FFTR1DInv subroutine uses only elements
    from 0th to floor(N/2)-th. It ignores imaginary part of F[0],  and in case
    N is even it ignores imaginary part of F[floor(N/2)] too.

    When you call this function using full arguments list - "FFTR1DInv(F,N,A)"
    - you can pass either either frequencies array with N elements or  reduced
    array with roughly N/2 elements - subroutine will  successfully  transform
    both.

    If you call this function using reduced arguments list -  "FFTR1DInv(F,A)"
    - you must pass FULL array with N elements (although higher  N/2 are still
    not used) because array size is used to automatically determine FFT length


      -- ALGLIB --
         Copyright 01.06.2009 by Bochkanov Sergey
    *************************************************************************/
    public static double[] FFT_Real_Inverse(Complex[] f, int n)
    {
        //
        // Special case: N=1, FFT is just identity transform.
        // After this block we assume that N is strictly greater than 1.
        //
        if (n == 1) return [f[0].Re];

        //
        // inverse real FFT is reduced to the inverse real FHT,
        // which is reduced to the forward real FHT,
        // which is reduced to the forward real FFT.
        //
        // Don't worry, it is really compact and efficient reduction :)
        //
        var h = new double[n];
        var result = new double[n];
        h[0] = f[0].Re;
        var n05 = (int)Floor(n / 2d);
        for (var i = 1; i < n05; i++)
        {
            var (re, im) = f[i];
            (h[i], h[n - i]) = (re - im, re + im);
        }

        if (n % 2 == 0) h[n05] = f[n05].Re;
        else
        {
            var (re, im) = f[n05];
            (h[n05], h[n05 + 1]) = (re - im, re + im);
        }

        var fh = FFT(h);

        for (var i = 0; i < n; i++)
        {
            var (re, im) = fh[i];
            result[i] = (re - im) / n;
        }
        return result;
    }


    /*************************************************************************
    Internal subroutine. Never call it directly!


      -- ALGLIB --
         Copyright 01.06.2009 by Bochkanov Sergey
    *************************************************************************/


    //private static void FFT_Real_InternalEven(ref double[] a, int n, ref double[] buf, FT_Base.ftplan plan)
    //{
    //    //
    //    // Special cases:
    //    // * N=2
    //    //
    //    // After this block we assume that N is strictly greater than 2
    //    //
    //    if (n == 2)
    //    {
    //        var x = a[0] + a[1];
    //        var y = a[0] - a[1];
    //        a[0] = x;
    //        a[1] = y;
    //        return;
    //    }

    //    //
    //    // even-size real FFT, use reduction to the Complex task
    //    //
    //    var n2 = n / 2;
    //    for (var i = 0; i < n; i++)
    //        buf[i] = a[i];

    //    FT_Base.FTBaseExecutePlan(ref buf, 0, plan);

    //    a[0] = buf[0] + buf[1];
    //    var pin = -Consts.pi2 / n;
    //    for (var i = 1; i < n2; i++)
    //    {
    //        var idx = 2 * (i % n2);
    //        var hn = new Complex(buf[idx + 0], buf[idx + 1]);
    //        idx = 2 * (n2 - i);
    //        var hmnc = new Complex(buf[idx + 0], -buf[idx + 1]);
    //        var v = new Complex(-Math.Sin(pin * i), Math.Cos(pin * i));
    //        v = hn + hmnc - v * (hn - hmnc);
    //        a[2 * i + 0] = .5 * v.Re;
    //        a[2 * i + 1] = .5 * v.Im;
    //    }
    //    a[1] = buf[0] - buf[1];
    //}



    /*************************************************************************
    Internal subroutine. Never call it directly!


      -- ALGLIB --
         Copyright 01.06.2009 by Bochkanov Sergey
    *************************************************************************/


    //private static void FFT_Real_Inverse_InternalEven(ref double[] a, int n, ref double[] buf, FT_Base.ftplan plan)
    //{
    //    //
    //    // Special cases:
    //    // * N=2
    //    //
    //    // After this block we assume that N is strictly greater than 2
    //    //
    //    if (n == 2)
    //    {
    //        var x = .5 * (a[0] + a[1]);
    //        var y = .5 * (a[0] - a[1]);
    //        a[0] = x;
    //        a[1] = y;
    //        return;
    //    }

    //    //
    //    // inverse real FFT is reduced to the inverse real FHT,
    //    // which is reduced to the forward real FHT,
    //    // which is reduced to the forward real FFT.
    //    //
    //    // Don't worry, it is really compact and efficient reduction :)
    //    //
    //    var n2 = n / 2;
    //    buf[0] = a[0];
    //    for (var i = 1; i < n2; i++)
    //    {
    //        var x = a[2 * i + 0];
    //        var y = a[2 * i + 1];
    //        buf[i] = x - y;
    //        buf[n - i] = x + y;
    //    }

    //    buf[n2] = a[1];
    //    FFT_Real_InternalEven(ref buf, n, ref a, plan);

    //    a[0] = buf[0] / n;
    //    var t = 1.0 / n;
    //    for (var i = 1; i < n2; i++)
    //    {
    //        var x = buf[2 * i + 0];
    //        var y = buf[2 * i + 1];
    //        a[i] = t * (x - y);
    //        a[n - i] = t * (x + y);
    //    }
    //    a[n2] = buf[1] / n;
    //}



    private static class FT_Base
    {
        public const int c_FT_BbasePlanEntrySize = 8;
        public const int c_FT_BaseCfftTask = 0;
        public const int ftbaserfhttask = 1;
        public const int FT_BaserFFT_Task = 2;
        public const int FFT_CooleyTukeyPlan = 0;
        public const int FFT_BluesteinPlan = 1;
        public const int FFT_CodeletPlan = 2;
        public const int FHT_CooleyTukeyPlan = 3;
        public const int FHT_CodeletPlan = 4;
        public const int FFT_RealCooleyTukeyPlan = 5;
        public const int FFT_EmptyPlan = 6;
        public const int FHT_N2Plan = 999;
        public const int FT_BaseUpdateTW = 4;
        public const int FT_BaseCodeletRecommended = 5;
        /*
                    public const double ftbaseinefficiencyfactor = 1.3;
        */
        public const int FT_BaseMaxSmoothFactor = 5;


        /*************************************************************************
        This subroutine generates FFT plan - a decomposition of a N-length FFT to
        the more simpler operations. Plan consists of the root entry and the child
        entries.

        Subroutine parameters:
            N               task size
        
        Output parameters:
            Plan            plan

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        public static void FTBaseGenerateComplexFFTPlan(int n, FT_Plan plan)
        {
            var plan_array_size = 1;
            plan.plan = new int[plan_array_size];

            var plan_size = 0;
            var pre_computed_size = 0;
            var stack_mem_size = 0;
            var stack_ptr = 0;
            var tmp_mem_size = 2 * n;
            FT_BaseGeneratePlanRec(
                n,
                c_FT_BaseCfftTask,
                plan,
                ref plan_size,
                ref pre_computed_size,
                ref plan_array_size,
                ref tmp_mem_size,
                ref stack_mem_size,
                stack_ptr);

            plan.StackBuffer = new double[Max(stack_mem_size, 1)];
            plan.TempBuffer = new double[Max(tmp_mem_size, 1)];
            plan.PreComputed = new double[Max(pre_computed_size, 1)];
            stack_ptr = 0;

            FT_BasePrecomputedPlanRec(plan, 0, stack_ptr);
        }


        /*************************************************************************
        Generates real FFT plan
        *************************************************************************/
        //public static void ftbasegeneraterealfftplan(int n, ftplan plan)
        //{
        //    var plan_array_size = 1;
        //    plan.plan = new int[plan_array_size];

        //    var plan_size = 0;
        //    var precomputed_size = 0;
        //    var stack_mem_size = 0;
        //    var stack_ptr = 0;
        //    var tmp_mem_size = 2 * n;
        //    FTBaseGeneratePlanRec(
        //        n,
        //        ft_base_rfft_task,
        //        plan,
        //        ref plan_size,
        //        ref precomputed_size,
        //        ref plan_array_size,
        //        ref tmp_mem_size,
        //        ref stack_mem_size,
        //        stack_ptr);

        //    plan.stackbuf = new double[Max(stack_mem_size, 1)];
        //    plan.tmpbuf = new double[Max(tmp_mem_size, 1)];
        //    plan.precomputed = new double[Max(precomputed_size, 1)];

        //    stack_ptr = 0;
        //    FTBasePrecomputePlanRec(plan, 0, stack_ptr);
        //}


        /*************************************************************************
        Generates real FHT plan
        *************************************************************************/
        //public static void ftbasegeneraterealfhtplan(int n, ftplan plan)
        //{
        //    var planarraysize = 1;
        //    plan.plan = new int[planarraysize];

        //    var plansize = 0;
        //    var precomputed_size = 0;
        //    var stack_mem_size = 0;
        //    var stack_ptr = 0;
        //    var tmp_mem_size = n;
        //    FTBaseGeneratePlanRec(n, ftbaserfhttask, plan, ref plansize, ref precomputed_size, ref planarraysize,
        //                          ref tmp_mem_size, ref stack_mem_size, stack_ptr);

        //    plan.stackbuf = new double[Math.Max(stack_mem_size, 1)];
        //    plan.tmpbuf = new double[Math.Max(tmp_mem_size, 1)];
        //    plan.precomputed = new double[Math.Max(precomputed_size, 1)];
        //    stack_ptr = 0;

        //    FTBasePrecomputePlanRec(plan, 0, stack_ptr);
        //}


        /*************************************************************************
        This subroutine executes FFT/FHT plan.

        If Plan is a:
        * complex FFT plan  -   sizeof(A)=2*N,
                                A contains interleaved real/imaginary values
        * real FFT plan     -   sizeof(A)=2*N,
                                A contains real values interleaved with zeros
        * real FHT plan     -   sizeof(A)=2*N,
                                A contains real values interleaved with zeros

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        public static void FTBaseExecutePlan(ref double[] a, int aoffset, FT_Plan plan)
        {
            const int stack_ptr = 0;
            FT_BaseExecutePlanRec(ref a, aoffset, plan, 0, stack_ptr);
        }


        /*************************************************************************
        Recurrent subroutine for the FTBaseExecutePlan

        Parameters:
            A           FFT'ed array
            AOffset     offset of the FFT'ed part (distance is measured in doubles)

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/

        private static void FT_BaseExecutePlanRec(ref double[] A, int Aoffset, FT_Plan plan, int EntryOffset, int StackPtr)
        {
            if (plan.plan[EntryOffset + 3] == FFT_EmptyPlan)
                return;

            int i;
            int n1;
            int n2;
            if (plan.plan[EntryOffset + 3] == FFT_CooleyTukeyPlan)
            {
                //
                // Cooley-Tukey plan
                // * transposition
                // * row-wise FFT
                // * twiddle factors:
                //   - TwBase is a basis twiddle factor for I=1, J=1
                //   - TwRow is a twiddle factor for a second element in a row (J=1)
                //   - Tw is a twiddle factor for a current element
                // * transposition again
                // * row-wise FFT again
                //
                n1 = plan.plan[EntryOffset + 1];
                n2 = plan.plan[EntryOffset + 2];

                InternalComplexLinTranspose(ref A, n1, n2, Aoffset, ref plan.TempBuffer);

                for (i = 0; i <= n2 - 1; i++)
                    FT_BaseExecutePlanRec(ref A, Aoffset + i * n1 * 2, plan, plan.plan[EntryOffset + 5], StackPtr);

                FFT_TW_Calc(ref A, Aoffset, n1, n2);
                InternalComplexLinTranspose(ref A, n2, n1, Aoffset, ref plan.TempBuffer);

                for (i = 0; i <= n1 - 1; i++)
                    FT_BaseExecutePlanRec(ref A, Aoffset + i * n2 * 2, plan, plan.plan[EntryOffset + 6], StackPtr);

                InternalComplexLinTranspose(ref A, n1, n2, Aoffset, ref plan.TempBuffer);
                return;
            }

            int offs;
            int offs2;
            double hk;
            double hnk;
            if (plan.plan[EntryOffset + 3] == FFT_RealCooleyTukeyPlan)
            {
                //
                // Cooley-Tukey plan
                // * transposition
                // * row-wise FFT
                // * twiddle factors:
                //   - TwBase is a basis twiddle factor for I=1, J=1
                //   - TwRow is a twiddle factor for a second element in a row (J=1)
                //   - Tw is a twiddle factor for a current element
                // * transposition again
                // * row-wise FFT again
                //
                n1 = plan.plan[EntryOffset + 1];
                n2 = plan.plan[EntryOffset + 2];

                InternalComplexLinTranspose(ref A, n2, n1, Aoffset, ref plan.TempBuffer);


                for (i = 0; i <= n1 / 2 - 1; i++)
                {
                    //
                    // pack two adjacent smaller real FFT's together,
                    // make one complex FFT,
                    // unpack result
                    //
                    offs = Aoffset + 2 * i * n2 * 2;
                    int k;
                    for (k = 0; k <= n2 - 1; k++)
                        A[offs + 2 * k + 1] = A[offs + 2 * n2 + 2 * k + 0];

                    FT_BaseExecutePlanRec(ref A, offs, plan, plan.plan[EntryOffset + 6], StackPtr);
                    plan.TempBuffer[0] = A[offs + 0];
                    plan.TempBuffer[1] = 0;
                    plan.TempBuffer[2 * n2 + 0] = A[offs + 1];
                    plan.TempBuffer[2 * n2 + 1] = 0;

                    for (k = 1; k <= n2 - 1; k++)
                    {
                        var offs1 = 2 * k;
                        offs2 = 2 * n2 + 2 * k;
                        hk = A[offs + 2 * k + 0];
                        hnk = A[offs + 2 * (n2 - k) + 0];
                        plan.TempBuffer[offs1 + 0] = 0.5 * (hk + hnk);
                        plan.TempBuffer[offs2 + 1] = -(0.5 * (hk - hnk));
                        hk = A[offs + 2 * k + 1];
                        hnk = A[offs + 2 * (n2 - k) + 1];
                        plan.TempBuffer[offs2 + 0] = .5 * (hk + hnk);
                        plan.TempBuffer[offs1 + 1] = .5 * (hk - hnk);
                    }

                    var i1_ = 0 - offs;
                    for (var i_ = offs; i_ <= offs + 2 * n2 * 2 - 1; i_++)
                        A[i_] = plan.TempBuffer[i_ + i1_];
                }
                if (n1 % 2 != 0)
                    FT_BaseExecutePlanRec(ref A, Aoffset + (n1 - 1) * n2 * 2, plan, plan.plan[EntryOffset + 6], StackPtr);

                FFT_TW_Calc(ref A, Aoffset, n2, n1);
                InternalComplexLinTranspose(ref A, n1, n2, Aoffset, ref plan.TempBuffer);

                for (i = 0; i <= n2 - 1; i++)
                    FT_BaseExecutePlanRec(ref A, Aoffset + i * n1 * 2, plan, plan.plan[EntryOffset + 5], StackPtr);

                InternalComplexLinTranspose(ref A, n2, n1, Aoffset, ref plan.TempBuffer);
                return;
            }

            int offsa;
            if (plan.plan[EntryOffset + 3] == FHT_CooleyTukeyPlan)
            {
                //
                // Cooley-Tukey FHT plan:
                // * transpose                    \
                // * smaller FHT's                |
                // * pre-process                  |
                // * multiply by twiddle factors  || corresponds to multiplication by H1
                // * post-process                 |
                // * transpose again              /
                // * multiply by H2 (smaller FHT's)
                // * final transposition
                //
                // For more details see Vitezslav Vesely, "Fast algorithms
                // of Fourier and Hartley transform and their implementation in MATLAB",
                // page 31.
                //
                n1 = plan.plan[EntryOffset + 1];
                n2 = plan.plan[EntryOffset + 2];
                //n = n1 * n2;

                InternalRealLinTranspose(ref A, n1, n2, Aoffset, ref plan.TempBuffer);

                for (i = 0; i <= n2 - 1; i++)
                    FT_BaseExecutePlanRec(ref A, Aoffset + i * n1, plan, plan.plan[EntryOffset + 5], StackPtr);
                int j;
                for (i = 0; i <= n2 - 1; i++)
                    for (j = 0; j <= n1 - 1; j++)
                    {
                        offsa = Aoffset + i * n1;
                        hk = A[offsa + j];
                        hnk = A[offsa + (n1 - j) % n1];
                        offs = 2 * (i * n1 + j);

                        plan.TempBuffer[offs + 0] = -.5 * (hnk - hk);
                        plan.TempBuffer[offs + 1] = .5 * (hk + hnk);
                    }

                FFT_TW_Calc(ref plan.TempBuffer, 0, n1, n2);

                for (j = 0; j <= n1 - 1; j++)
                    A[Aoffset + j] = plan.TempBuffer[2 * j + 0] + plan.TempBuffer[2 * j + 1];

                if (n2 % 2 == 0)
                {
                    offs = 2 * (n2 / 2) * n1;
                    offsa = Aoffset + n2 / 2 * n1;

                    for (j = 0; j <= n1 - 1; j++)
                        A[offsa + j] = plan.TempBuffer[offs + 2 * j + 0] + plan.TempBuffer[offs + 2 * j + 1];
                }

                for (i = 1; i <= (n2 + 1) / 2 - 1; i++)
                {
                    offs = 2 * i * n1;
                    offs2 = 2 * (n2 - i) * n1;
                    offsa = Aoffset + i * n1;

                    for (j = 0; j <= n1 - 1; j++)
                        A[offsa + j] = plan.TempBuffer[offs + 2 * j + 1] + plan.TempBuffer[offs2 + 2 * j + 0];

                    offsa = Aoffset + (n2 - i) * n1;

                    for (j = 0; j <= n1 - 1; j++)
                        A[offsa + j] = plan.TempBuffer[offs + 2 * j + 0] + plan.TempBuffer[offs2 + 2 * j + 1];
                }

                InternalRealLinTranspose(ref A, n2, n1, Aoffset, ref plan.TempBuffer);

                for (i = 0; i <= n1 - 1; i++)
                    FT_BaseExecutePlanRec(ref A, Aoffset + i * n2, plan, plan.plan[EntryOffset + 6], StackPtr);

                InternalRealLinTranspose(ref A, n1, n2, Aoffset, ref plan.TempBuffer);
                return;
            }

            int n;
            if (plan.plan[EntryOffset + 3] == FHT_N2Plan)
            {
                //
                // Cooley-Tukey FHT plan
                //
                n1 = plan.plan[EntryOffset + 1];
                n2 = plan.plan[EntryOffset + 2];
                n = n1 * n2;

                RefFHT(ref A, n, Aoffset);
                return;
            }

            double s1x;
            double s2x;
            double s3y;
            double s4x;
            double s5y;

            double c1;
            double c2;
            double c3;
            double c4;
            double c5;

            double m1x;
            double m2x;
            double m2y;
            double m3y;

            double t1x;
            double t2x;
            double t3x;
            double t4x;
            double t5x;

            double a0x;
            double a1x;
            double a2x;
            double a3x;

            double v0;
            if (plan.plan[EntryOffset + 3] == FFT_CodeletPlan)
            {
                n1 = plan.plan[EntryOffset + 1];
                n2 = plan.plan[EntryOffset + 2];
                n = n1 * n2;

                double a0y;
                double a1y;
                if (n == 2)
                {
                    a0x = A[Aoffset + 0];
                    a0y = A[Aoffset + 1];
                    a1x = A[Aoffset + 2];
                    a1y = A[Aoffset + 3];

                    v0 = a0x + a1x;
                    var v1 = a0y + a1y;
                    var v2 = a0x - a1x;
                    var v3 = a0y - a1y;

                    A[Aoffset + 0] = v0;
                    A[Aoffset + 1] = v1;
                    A[Aoffset + 2] = v2;
                    A[Aoffset + 3] = v3;
                    return;
                }

                double s1y;
                double t1y;
                double a2y;


                double m1y;
                if (n == 3)
                {
                    offs = plan.plan[EntryOffset + 7];

                    c1 = plan.PreComputed[offs + 0];
                    c2 = plan.PreComputed[offs + 1];

                    a0x = A[Aoffset + 0];
                    a0y = A[Aoffset + 1];
                    a1x = A[Aoffset + 2];
                    a1y = A[Aoffset + 3];
                    a2x = A[Aoffset + 4];
                    a2y = A[Aoffset + 5];

                    t1x = a1x + a2x;
                    t1y = a1y + a2y;

                    a0x += t1x;
                    a0y += t1y;

                    m1x = c1 * t1x;
                    m1y = c1 * t1y;
                    m2x = c2 * (a1y - a2y);
                    m2y = c2 * (a2x - a1x);

                    s1x = a0x + m1x;
                    s1y = a0y + m1y;

                    a1x = s1x + m2x;
                    a1y = s1y + m2y;
                    a2x = s1x - m2x;
                    a2y = s1y - m2y;

                    A[Aoffset + 0] = a0x;
                    A[Aoffset + 1] = a0y;
                    A[Aoffset + 2] = a1x;
                    A[Aoffset + 3] = a1y;
                    A[Aoffset + 4] = a2x;
                    A[Aoffset + 5] = a2y;
                    return;
                }

                double t2y;
                double m3x;
                if (n == 4)
                {
                    a0x = A[Aoffset + 0];
                    a0y = A[Aoffset + 1];
                    a1x = A[Aoffset + 2];
                    a1y = A[Aoffset + 3];
                    a2x = A[Aoffset + 4];
                    a2y = A[Aoffset + 5];
                    a3x = A[Aoffset + 6];
                    var a3y = A[Aoffset + 7];

                    t1x = a0x + a2x;
                    t1y = a0y + a2y;
                    t2x = a1x + a3x;
                    t2y = a1y + a3y;

                    m2x = a0x - a2x;
                    m2y = a0y - a2y;
                    m3x = a1y - a3y;
                    m3y = a3x - a1x;

                    A[Aoffset + 0] = t1x + t2x;
                    A[Aoffset + 1] = t1y + t2y;
                    A[Aoffset + 4] = t1x - t2x;
                    A[Aoffset + 5] = t1y - t2y;
                    A[Aoffset + 2] = m2x + m3x;
                    A[Aoffset + 3] = m2y + m3y;
                    A[Aoffset + 6] = m2x - m3x;
                    A[Aoffset + 7] = m2y - m3y;
                    return;
                }

                if (n == 5)
                {
                    offs = plan.plan[EntryOffset + 7];

                    c1 = plan.PreComputed[offs + 0];
                    c2 = plan.PreComputed[offs + 1];
                    c3 = plan.PreComputed[offs + 2];
                    c4 = plan.PreComputed[offs + 3];
                    c5 = plan.PreComputed[offs + 4];

                    t1x = A[Aoffset + 2] + A[Aoffset + 8];
                    t1y = A[Aoffset + 3] + A[Aoffset + 9];
                    t2x = A[Aoffset + 4] + A[Aoffset + 6];
                    t2y = A[Aoffset + 5] + A[Aoffset + 7];
                    t3x = A[Aoffset + 2] - A[Aoffset + 8];
                    var t3y = A[Aoffset + 3] - A[Aoffset + 9];
                    t4x = A[Aoffset + 6] - A[Aoffset + 4];
                    var t4y = A[Aoffset + 7] - A[Aoffset + 5];

                    t5x = t1x + t2x;
                    var t5y = t1y + t2y;

                    A[Aoffset + 0] += t5x;
                    A[Aoffset + 1] += t5y;

                    m1x = c1 * t5x;
                    m1y = c1 * t5y;
                    m2x = c2 * (t1x - t2x);
                    m2y = c2 * (t1y - t2y);
                    m3x = -(c3 * (t3y + t4y));
                    m3y = c3 * (t3x + t4x);
                    var m4x = -(c4 * t4y);
                    var m4y = c4 * t4x;
                    var m5x = -(c5 * t3y);
                    var m5y = c5 * t3x;

                    var s3x = m3x - m4x;
                    s3y = m3y - m4y;

                    var s5x = m3x + m5x;
                    s5y = m3y + m5y;

                    s1x = A[Aoffset + 0] + m1x;
                    s1y = A[Aoffset + 1] + m1y;
                    s2x = s1x + m2x;
                    var s2y = s1y + m2y;
                    s4x = s1x - m2x;
                    var s4y = s1y - m2y;

                    A[Aoffset + 2] = s2x + s3x;
                    A[Aoffset + 3] = s2y + s3y;
                    A[Aoffset + 4] = s4x + s5x;
                    A[Aoffset + 5] = s4y + s5y;
                    A[Aoffset + 6] = s4x - s5x;
                    A[Aoffset + 7] = s4y - s5y;
                    A[Aoffset + 8] = s2x - s3x;
                    A[Aoffset + 9] = s2y - s3y;
                    return;
                }
            }

            if (plan.plan[EntryOffset + 3] == FHT_CodeletPlan)
            {
                n1 = plan.plan[EntryOffset + 1];
                n2 = plan.plan[EntryOffset + 2];
                n = n1 * n2;

                if (n == 2)
                {
                    a0x = A[Aoffset + 0];
                    a1x = A[Aoffset + 1];

                    A[Aoffset + 0] = a0x + a1x;
                    A[Aoffset + 1] = a0x - a1x;
                    return;
                }

                if (n == 3)
                {
                    offs = plan.plan[EntryOffset + 7];

                    c1 = plan.PreComputed[offs + 0];
                    c2 = plan.PreComputed[offs + 1];

                    a0x = A[Aoffset + 0];
                    a1x = A[Aoffset + 1];
                    a2x = A[Aoffset + 2];

                    t1x = a1x + a2x;
                    a0x += t1x;

                    m1x = c1 * t1x;
                    m2y = c2 * (a2x - a1x);

                    s1x = a0x + m1x;

                    A[Aoffset + 0] = a0x;
                    A[Aoffset + 1] = s1x - m2y;
                    A[Aoffset + 2] = s1x + m2y;
                    return;
                }

                if (n == 4)
                {
                    a0x = A[Aoffset + 0];
                    a1x = A[Aoffset + 1];
                    a2x = A[Aoffset + 2];
                    a3x = A[Aoffset + 3];

                    t1x = a0x + a2x;
                    t2x = a1x + a3x;

                    m2x = a0x - a2x;
                    m3y = a3x - a1x;

                    A[Aoffset + 0] = t1x + t2x;
                    A[Aoffset + 1] = m2x - m3y;
                    A[Aoffset + 2] = t1x - t2x;
                    A[Aoffset + 3] = m2x + m3y;
                    return;
                }

                if (n == 5)
                {
                    offs = plan.plan[EntryOffset + 7];

                    c1 = plan.PreComputed[offs + 0];
                    c2 = plan.PreComputed[offs + 1];
                    c3 = plan.PreComputed[offs + 2];
                    c4 = plan.PreComputed[offs + 3];
                    c5 = plan.PreComputed[offs + 4];

                    t1x = A[Aoffset + 1] + A[Aoffset + 4];
                    t2x = A[Aoffset + 2] + A[Aoffset + 3];
                    t3x = A[Aoffset + 1] - A[Aoffset + 4];
                    t4x = A[Aoffset + 3] - A[Aoffset + 2];

                    t5x = t1x + t2x;
                    v0 = A[Aoffset + 0] + t5x;
                    A[Aoffset + 0] = v0;

                    m2x = c2 * (t1x - t2x);
                    m3y = c3 * (t3x + t4x);

                    s3y = m3y - c4 * t4x;
                    s5y = m3y + c5 * t3x;
                    s1x = v0 + c1 * t5x;

                    s2x = s1x + m2x;
                    s4x = s1x - m2x;

                    A[Aoffset + 1] = s2x - s3y;
                    A[Aoffset + 2] = s4x - s5y;
                    A[Aoffset + 3] = s4x + s5y;
                    A[Aoffset + 4] = s2x + s3y;
                    return;
                }
            }

            if (plan.plan[EntryOffset + 3] != FFT_BluesteinPlan) return;
            //
            // Bluestein plan:
            // 1. multiply by precomputed coefficients
            // 2. make convolution: forward FFT, multiplication by precomputed FFT
            //    and backward FFT. backward FFT is represented as
            //
            //        invfft(x) = fft(x')'/M
            //
            //    for performance reasons reduction of inverse FFT to
            //    forward FFT is merged with multiplication of FFT components
            //    and last stage of Bluestein's transformation.
            // 3. post-multiplication by Bluestein factors
            //
            n = plan.plan[EntryOffset + 1];
            var m = plan.plan[EntryOffset + 4];
            offs = plan.plan[EntryOffset + 7];

            for (i = StackPtr + 2 * n; i <= StackPtr + 2 * m - 1; i++)
                plan.StackBuffer[i] = 0;

            var offsp = offs + 2 * m;
            offsa = Aoffset;
            var offsb = StackPtr;

            double x;
            double y;
            double bx;
            double by;
            for (i = 0; i <= n - 1; i++)
            {
                bx = plan.PreComputed[offsp + 0];
                by = plan.PreComputed[offsp + 1];

                x = A[offsa + 0];
                y = A[offsa + 1];

                plan.StackBuffer[offsb + 0] = x * bx - y * -by;
                plan.StackBuffer[offsb + 1] = x * -by + y * bx;

                offsp += 2;
                offsa += 2;
                offsb += 2;
            }

            FT_BaseExecutePlanRec(ref plan.StackBuffer, StackPtr, plan, plan.plan[EntryOffset + 5], StackPtr + 2 * 2 * m);

            offsb = StackPtr;
            offsp = offs;
            for (i = 0; i <= m - 1; i++)
            {
                x = plan.StackBuffer[offsb + 0];
                y = plan.StackBuffer[offsb + 1];

                bx = plan.PreComputed[offsp + 0];
                by = plan.PreComputed[offsp + 1];

                plan.StackBuffer[offsb + 0] = x * bx - y * by;
                plan.StackBuffer[offsb + 1] = -(x * by + y * bx);

                offsb += 2;
                offsp += 2;
            }

            FT_BaseExecutePlanRec(ref plan.StackBuffer, StackPtr, plan, plan.plan[EntryOffset + 5], StackPtr + 2 * 2 * m);

            offsb = StackPtr;
            offsp = offs + 2 * m;
            offsa = Aoffset;

            for (i = 0; i <= n - 1; i++)
            {
                x = plan.StackBuffer[offsb + 0] / m;
                y = -(plan.StackBuffer[offsb + 1] / m);

                bx = plan.PreComputed[offsp + 0];
                by = plan.PreComputed[offsp + 1];

                A[offsa + 0] = x * bx - y * -by;
                A[offsa + 1] = x * -by + y * bx;

                offsp += 2;
                offsa += 2;
                offsb += 2;
            }
        }

        /*************************************************************************
        Returns good factorization N=N1*N2.

        Usually N1<=N2 (but not always - small N's may be exception).
        if N1<>1 then N2<>1.

        Factorization is chosen depending on task type and codelets we have.

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/

        private static void FT_BaseFactorize(int n, out int n1, out int n2)
        {
            int j;

            n1 = 0;
            n2 = 0;

            //
            // try to find good codelet
            //
            if (n1 * n2 != n)
                for (j = FT_BaseCodeletRecommended; j >= 2; j--)
                    if (n % j == 0)
                    {
                        n1 = j;
                        n2 = n / j;
                        break;
                    }

            //
            // try to factorize N
            //
            if (n1 * n2 != n)
                for (j = FT_BaseCodeletRecommended + 1; j <= n - 1; j++)
                    if (n % j == 0)
                    {
                        n1 = j;
                        n2 = n / j;
                        break;
                    }

            //
            // looks like N is prime :(
            //
            if (n1 * n2 != n)
            {
                n1 = 1;
                n2 = n;
            }

            //
            // normalize
            //
            if (n2 != 1 || n1 == 1) return;
            n2 = n1;
            n1 = 1;
        }


        /*************************************************************************
        Is number smooth?

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        /*
                    public static bool ftbaseissmooth(int n)
                    {

                        for(var i = 2; i <= ftbasemaxsmoothfactor; i++)
                            while(n % i == 0)
                                n /= i;
                        return n == 1;
                    }
        */


        /*************************************************************************
        Returns smallest smooth (divisible only by 2, 3, 5) number that is greater
        than or equal to max(N,2)

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/

        private static int FT_BaseFindSmooth(int n)
        {
            var best = 2;
            while (best < n) best = 2 * best;
            FT_BaseFindSmoothRec(n, 1, 2, ref best);
            return best;
        }


        /*************************************************************************
        Returns  smallest  smooth  (divisible only by 2, 3, 5) even number that is
        greater than or equal to max(N,2)

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        /*
                    public static int ftbasefindsmootheven(int n)
                    {
                        var best = 2;
                        while(best < n) best = 2 * best;
                        ftbasefindsmoothrec(n, 2, 2, ref best);
                        return best;
                    }
        */


        /*************************************************************************
        Returns estimate of FLOP count for the FFT.

        It is only an estimate based on operations count for the PERFECT FFT
        and relative inefficiency of the algorithm actually used.

        N should be power of 2, estimates are badly wrong for non-power-of-2 N's.

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        /*
                    public static double ftbasegetflopestimate(int n)
                    {
                        return ftbaseinefficiencyfactor * (4 * n * Math.Log(n) / Math.Log(2) - 6 * n + 8);
                    }
        */


        /*************************************************************************
        Recurrent subroutine for the FFTGeneratePlan:

        PARAMETERS:
            N                   plan size
            IsReal              whether input is real or not.
                                subroutine MUST NOT ignore this flag because real
                                inputs comes with non-initialized imaginary parts,
                                so ignoring this flag will result in corrupted output
            HalfOut             whether full output or only half of it from 0 to
                                floor(N/2) is needed. This flag may be ignored if
                                doing so will simplify calculations
            Plan                plan array
            PlanSize            size of used part (in integers)
            PrecomputedSize     size of precomputed array allocated yet
            PlanArraySize       plan array size (actual)
            TmpMemSize          temporary memory required size
            BluesteinMemSize    temporary memory required size

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        private static void FT_BaseGeneratePlanRec
        (
            int n,
            int TaskType,
            FT_Plan plan,
            ref int PlanSize,
            ref int PreComputedSize,
            ref int PlanArraySize,
            ref int TempMemorySize,
            ref int StackMemorySize,
            int StackPtr
        )
        {
            //
            // prepare
            //
            if (PlanSize + c_FT_BbasePlanEntrySize > PlanArraySize)
                FFT_AarrayResize(ref plan.plan, ref PlanArraySize, 8 * PlanArraySize);
            var entryoffset = PlanSize;
            const int esize = c_FT_BbasePlanEntrySize;
            PlanSize += esize;

            //
            // if N=1, generate empty plan and exit
            //
            if (n == 1)
            {
                plan.plan[entryoffset + 0] = esize;
                plan.plan[entryoffset + 1] = -1;
                plan.plan[entryoffset + 2] = -1;
                plan.plan[entryoffset + 3] = FFT_EmptyPlan;
                plan.plan[entryoffset + 4] = -1;
                plan.plan[entryoffset + 5] = -1;
                plan.plan[entryoffset + 6] = -1;
                plan.plan[entryoffset + 7] = -1;
                return;
            }

            //
            // generate plans
            //
            FT_BaseFactorize(n, out var n1, out var n2);
            if (TaskType is c_FT_BaseCfftTask or FT_BaserFFT_Task)
            {
                //
                // complex FFT plans
                //
                if (n1 != 1)
                {
                    //
                    // Cooley-Tukey plan (real or complex)
                    //
                    // Note that child plans are COMPLEX
                    // (whether plan itself is complex or not).
                    //
                    TempMemorySize = Max(TempMemorySize, 2 * n1 * n2);
                    plan.plan[entryoffset + 0] = esize;
                    plan.plan[entryoffset + 1] = n1;
                    plan.plan[entryoffset + 2] = n2;

                    plan.plan[entryoffset + 3] = TaskType == c_FT_BaseCfftTask
                        ? FFT_CooleyTukeyPlan
                        : FFT_RealCooleyTukeyPlan;

                    plan.plan[entryoffset + 4] = 0;
                    plan.plan[entryoffset + 5] = PlanSize;
                    FT_BaseGeneratePlanRec(n1, c_FT_BaseCfftTask, plan, ref PlanSize, ref PreComputedSize, ref PlanArraySize,
                        ref TempMemorySize, ref StackMemorySize, StackPtr);
                    plan.plan[entryoffset + 6] = PlanSize;
                    FT_BaseGeneratePlanRec(n2, c_FT_BaseCfftTask, plan, ref PlanSize, ref PreComputedSize, ref PlanArraySize,
                        ref TempMemorySize, ref StackMemorySize, StackPtr);
                    plan.plan[entryoffset + 7] = -1;
                    return;
                }

                if (n is 2 or 3 or 4 or 5)
                {
                    //
                    // hard-coded plan
                    //
                    plan.plan[entryoffset + 0] = esize;
                    plan.plan[entryoffset + 1] = n1;
                    plan.plan[entryoffset + 2] = n2;
                    plan.plan[entryoffset + 3] = FFT_CodeletPlan;
                    plan.plan[entryoffset + 4] = 0;
                    plan.plan[entryoffset + 5] = -1;
                    plan.plan[entryoffset + 6] = -1;
                    plan.plan[entryoffset + 7] = PreComputedSize;

                    if (n == 3) PreComputedSize += 2;

                    if (n == 5) PreComputedSize += 5;
                    return;
                }
                //
                // Bluestein's plan
                //
                // Select such M that M>=2*N-1, M is composite, and M's
                // factors are 2, 3, 5
                //
                var k = 2 * n2 - 1;
                var m = FT_BaseFindSmooth(k);
                TempMemorySize = Max(TempMemorySize, 2 * m);
                plan.plan[entryoffset + 0] = esize;
                plan.plan[entryoffset + 1] = n2;
                plan.plan[entryoffset + 2] = -1;
                plan.plan[entryoffset + 3] = FFT_BluesteinPlan;
                plan.plan[entryoffset + 4] = m;
                plan.plan[entryoffset + 5] = PlanSize;

                StackPtr += 2 * 2 * m;
                StackMemorySize = Max(StackMemorySize, StackPtr);

                FT_BaseGeneratePlanRec(m, c_FT_BaseCfftTask, plan, ref PlanSize, ref PreComputedSize,
                    ref PlanArraySize, ref TempMemorySize, ref StackMemorySize, StackPtr);

                //stack_ptr = stack_ptr - 2 * 2 * m;
                plan.plan[entryoffset + 6] = -1;
                plan.plan[entryoffset + 7] = PreComputedSize;
                PreComputedSize = PreComputedSize + 2 * m + 2 * n;
                return;
            }
            if (TaskType != ftbaserfhttask) return;
            //
            // real FHT plans
            //
            if (n1 != 1)
            {
                //
                // Cooley-Tukey plan
                //
                //
                TempMemorySize = Max(TempMemorySize, 2 * n1 * n2);
                plan.plan[entryoffset + 0] = esize;
                plan.plan[entryoffset + 1] = n1;
                plan.plan[entryoffset + 2] = n2;
                plan.plan[entryoffset + 3] = FHT_CooleyTukeyPlan;
                plan.plan[entryoffset + 4] = 0;
                plan.plan[entryoffset + 5] = PlanSize;

                FT_BaseGeneratePlanRec(n1, TaskType, plan, ref PlanSize, ref PreComputedSize, ref PlanArraySize,
                    ref TempMemorySize, ref StackMemorySize, StackPtr);

                plan.plan[entryoffset + 6] = PlanSize;

                FT_BaseGeneratePlanRec(n2, TaskType, plan, ref PlanSize, ref PreComputedSize, ref PlanArraySize,
                    ref TempMemorySize, ref StackMemorySize, StackPtr);

                plan.plan[entryoffset + 7] = -1;
                return;
            }
            //
            // N2 plan
            //
            plan.plan[entryoffset + 0] = esize;
            plan.plan[entryoffset + 1] = n1;
            plan.plan[entryoffset + 2] = n2;
            plan.plan[entryoffset + 3] = FHT_N2Plan;
            plan.plan[entryoffset + 4] = 0;
            plan.plan[entryoffset + 5] = -1;
            plan.plan[entryoffset + 6] = -1;
            plan.plan[entryoffset + 7] = -1;

            if (n is 2 or 3 or 4 or 5)
            {
                //
                // hard-coded plan
                //
                plan.plan[entryoffset + 0] = esize;
                plan.plan[entryoffset + 1] = n1;
                plan.plan[entryoffset + 2] = n2;
                plan.plan[entryoffset + 3] = FHT_CodeletPlan;
                plan.plan[entryoffset + 4] = 0;
                plan.plan[entryoffset + 5] = -1;
                plan.plan[entryoffset + 6] = -1;
                plan.plan[entryoffset + 7] = PreComputedSize;

                if (n == 3) PreComputedSize += 2;

                if (n == 5) PreComputedSize += 5;
            }
        }


        /*************************************************************************
        Recurrent subroutine for precomputing FFT plans

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        private static void FT_BasePrecomputedPlanRec(FT_Plan plan, int entryoffset, int stack_ptr)
        {

            if (plan.plan[entryoffset + 3] == FFT_CooleyTukeyPlan ||
                plan.plan[entryoffset + 3] == FFT_RealCooleyTukeyPlan ||
                plan.plan[entryoffset + 3] == FHT_CooleyTukeyPlan)
            {
                FT_BasePrecomputedPlanRec(plan, plan.plan[entryoffset + 5], stack_ptr);
                FT_BasePrecomputedPlanRec(plan, plan.plan[entryoffset + 6], stack_ptr);
                return;
            }

            int n;
            int offs;
            if (plan.plan[entryoffset + 3] == FFT_CodeletPlan || plan.plan[entryoffset + 3] == FHT_CodeletPlan)
            {
                var n1 = plan.plan[entryoffset + 1];
                var n2 = plan.plan[entryoffset + 2];
                n = n1 * n2;

                const double pi23 = Consts.pi2 / 3;
                if (n == 3)
                {
                    offs = plan.plan[entryoffset + 7];
                    plan.PreComputed[offs + 0] = Cos(pi23) - 1;
                    plan.PreComputed[offs + 1] = Sin(pi23);
                    return;
                }

                if (n == 5)
                {
                    offs = plan.plan[entryoffset + 7];
                    const double pi25 = Consts.pi2 / 5;
                    const double pi45 = pi25 * 2;
                    plan.PreComputed[offs + 0] = (Cos(pi25) + Cos(pi45)) / 2 - 1;
                    plan.PreComputed[offs + 1] = (Cos(pi25) - Cos(pi45)) / 2;
                    plan.PreComputed[offs + 2] = -Sin(pi25);
                    plan.PreComputed[offs + 3] = -(Sin(pi25) + Sin(pi45));
                    plan.PreComputed[offs + 4] = Sin(pi25) - Sin(pi45);
                    return;
                }
            }

            if (plan.plan[entryoffset + 3] != FFT_BluesteinPlan) return;

            FT_BasePrecomputedPlanRec(plan, plan.plan[entryoffset + 5], stack_ptr);

            n = plan.plan[entryoffset + 1];
            var m = plan.plan[entryoffset + 4];
            offs = plan.plan[entryoffset + 7];

            int i;
            for (i = 0; i < 2 * m; i++)
                plan.PreComputed[offs + i] = 0;

            var pin = PI / n;
            for (i = 0; i < n; i++)
            {
                var bx = Cos(pin * i * i);
                var by = Sin(pin * i * i);
                plan.PreComputed[offs + 2 * i + 0] = bx;
                plan.PreComputed[offs + 2 * i + 1] = by;
                plan.PreComputed[offs + 2 * m + 2 * i + 0] = bx;
                plan.PreComputed[offs + 2 * m + 2 * i + 1] = by;

                if (i <= 0) continue;
                plan.PreComputed[offs + 2 * (m - i) + 0] = bx;
                plan.PreComputed[offs + 2 * (m - i) + 1] = by;
            }

            FT_BaseExecutePlanRec(ref plan.PreComputed, offs, plan, plan.plan[entryoffset + 5], stack_ptr);
        }


        /*************************************************************************
        Twiddle factors calculation

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        private static void FFT_TW_Calc(ref double[] a, int aoffset, int n1, int n2)
        {
            var n = n1 * n2;
            var pin = -Consts.pi2 / n;
            var twbasexm1 = Sin(.5 * pin);
            twbasexm1 *= -2 * twbasexm1;
            //twbasexm1 = -(2 * math.sqr(Math.Sin(0.5 * v)));
            var twbasey = Sin(pin);
            double twrowxm1 = 0;
            double twrowy = 0;
            for (var i = 0; i <= n2 - 1; i++)
            {
                double twxm1 = 0;
                double twy = 0;
                double tmpy;
                double tmpx;
                int j;
                for (j = 0; j <= n1 - 1; j++)
                {
                    var idx = i * n1 + j;
                    var offs = aoffset + 2 * idx;
                    var x = a[offs + 0];
                    var y = a[offs + 1];
                    tmpx = x * twxm1 - y * twy;
                    tmpy = x * twy + y * twxm1;
                    a[offs + 0] = x + tmpx;
                    a[offs + 1] = y + tmpy;

                    //
                    // update Tw: Tw(new) = Tw(old)*TwRow
                    //
                    if (j >= n1 - 1) continue;
                    if (j % FT_BaseUpdateTW == 0)
                    {
                        pin = -(Consts.pi2 * i * (j + 1) / n);
                        twxm1 = Sin(.5 * pin);
                        twxm1 *= -2 * twxm1;
                        //twxm1 = -(2 * math.sqr(Math.Sin(0.5 * v)));
                        twy = Sin(pin);
                    }
                    else
                    {
                        tmpx = twrowxm1 + twxm1 * twrowxm1 - twy * twrowy;
                        tmpy = twrowy + twxm1 * twrowy + twy * twrowxm1;
                        twxm1 += tmpx;
                        twy += tmpy;
                    }
                }

                //
                // update TwRow: TwRow(new) = TwRow(old)*TwBase
                //
                if (i >= n2 - 1) continue;
                if (j % FT_BaseUpdateTW == 0)
                {
                    pin = -(Consts.pi2 * (i + 1) / n);
                    twrowxm1 = Sin(.5 * pin);
                    twrowxm1 *= -2 * twrowxm1;
                    //twrowxm1 = -(2 * math.sqr(Math.Sin(0.5 * v)));
                    twrowy = Sin(pin);
                }
                else
                {
                    tmpx = twbasexm1 + twrowxm1 * twbasexm1 - twrowy * twbasey;
                    tmpy = twbasey + twrowxm1 * twbasey + twrowy * twbasexm1;
                    twrowxm1 += tmpx;
                    twrowy += tmpy;
                }
            }
        }


        /*************************************************************************
        Linear transpose: transpose complex matrix stored in 1-dimensional array

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        private static void InternalComplexLinTranspose(ref double[] a, int m, int n, int astart, ref double[] buf)
        {
            FFTicltRec(ref a, astart, n, ref buf, 0, m, m, n);
            var i1 = 0 - astart;
            for (var i = astart; i <= astart + 2 * m * n - 1; i++)
                a[i] = buf[i + i1];
        }


        /*************************************************************************
        Linear transpose: transpose real matrix stored in 1-dimensional array

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        private static void InternalRealLinTranspose(ref double[] a, int m, int n, int astart, ref double[] buf)
        {

            FFTirltRec(ref a, astart, n, ref buf, 0, m, m, n);
            var i1 = 0 - astart;
            for (var i = astart; i <= astart + m * n - 1; i++)
                a[i] = buf[i + i1];
        }


        /*************************************************************************
        Recurrent subroutine for a InternalComplexLinTranspose

        Write A^T to B, where:
        * A is m*n complex matrix stored in array A as pairs of real/image values,
          beginning from AStart position, with AStride stride
        * B is n*m complex matrix stored in array B as pairs of real/image values,
          beginning from BStart position, with BStride stride
        stride is measured in complex numbers, i.e. in real/image pairs.

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        private static void FFTicltRec(
            ref double[] a,
            int astart,
            int astride,
            ref double[] b,
            int bstart,
            int bstride,
            int m,
            int n)
        {
            if (m == 0 || n == 0) return;
            if (Max(m, n) <= 8)
            {
                var m2 = 2 * bstride;
                for (var i = 0; i <= m - 1; i++)
                {
                    var idx1 = bstart + 2 * i;
                    var idx2 = astart + 2 * i * astride;
                    for (var j = 0; j <= n - 1; j++)
                    {
                        b[idx1 + 0] = a[idx2 + 0];
                        b[idx1 + 1] = a[idx2 + 1];
                        idx1 += m2;
                        idx2 += 2;
                    }
                }
                return;
            }
            if (n > m)
            {
                //
                // New partition:
                //
                // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
                //                                  ( B2 )
                //
                var n1 = n / 2;
                if (n - n1 >= 8 && n1 % 8 != 0)
                    n1 += 8 - n1 % 8;
                FFTicltRec(ref a, astart, astride, ref b, bstart, bstride, m, n1);
                FFTicltRec(ref a, astart + 2 * n1, astride, ref b, bstart + 2 * n1 * bstride, bstride, m, n - n1);
            }
            else
            {
                //
                // New partition:
                //
                // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
                //                     ( A2 )
                //
                var m1 = m / 2;
                if (m - m1 >= 8 && m1 % 8 != 0)
                    m1 += 8 - m1 % 8;
                FFTicltRec(ref a, astart, astride, ref b, bstart, bstride, m1, n);
                FFTicltRec(ref a, astart + 2 * m1 * astride, astride, ref b, bstart + 2 * m1, bstride, m - m1, n);
            }
        }


        /*************************************************************************
        Recurrent subroutine for a InternalRealLinTranspose


          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        private static void FFTirltRec(ref double[] a, int astart, int astride, ref double[] b, int bstart,
            int bstride, int m, int n)
        {
            if (m == 0 || n == 0) return;

            if (Max(m, n) <= 8)
            {
                for (var i = 0; i <= m - 1; i++)
                {
                    var idx1 = bstart + i;
                    var idx2 = astart + i * astride;
                    for (var j = 0; j <= n - 1; j++)
                    {
                        b[idx1] = a[idx2];
                        idx1 += bstride;
                        idx2 += 1;
                    }
                }
                return;
            }
            if (n > m)
            {
                //
                // New partition:
                //
                // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
                //                                  ( B2 )
                //
                var n1 = n / 2;
                if (n - n1 >= 8 && n1 % 8 != 0)
                    n1 += 8 - n1 % 8;
                FFTirltRec(ref a, astart, astride, ref b, bstart, bstride, m, n1);
                FFTirltRec(ref a, astart + n1, astride, ref b, bstart + n1 * bstride, bstride, m, n - n1);
            }
            else
            {
                //
                // New partition:
                //
                // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
                //                     ( A2 )
                //
                var m1 = m / 2;
                if (m - m1 >= 8 && m1 % 8 != 0)
                    m1 += 8 - m1 % 8;
                FFTirltRec(ref a, astart, astride, ref b, bstart, bstride, m1, n);
                FFTirltRec(ref a, astart + m1 * astride, astride, ref b, bstart + m1, bstride, m - m1, n);
            }
        }


        /*************************************************************************
        recurrent subroutine for FFTFindSmoothRec

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        private static void FT_BaseFindSmoothRec(int n, int seed, int leastfactor, ref int best)
        {
            if (seed >= n)
            {
                best = Min(best, seed);
                return;
            }
            if (leastfactor <= 2)
                FT_BaseFindSmoothRec(n, seed * 2, 2, ref best);
            if (leastfactor <= 3)
                FT_BaseFindSmoothRec(n, seed * 3, 3, ref best);
            if (leastfactor <= 5)
                FT_BaseFindSmoothRec(n, seed * 5, 5, ref best);
        }


        /*************************************************************************
        Internal subroutine: array resize

          -- ALGLIB --
             Copyright 01.05.2009 by Bochkanov Sergey
        *************************************************************************/
        private static void FFT_AarrayResize(ref int[] a, ref int asize, int newasize)
        {
            var tmp = new int[asize];
            for (var i = 0; i <= asize - 1; i++)
                tmp[i] = a[i];

            a = new int[newasize];
            for (var i = 0; i <= asize - 1; i++)
                a[i] = tmp[i];

            asize = newasize;
        }


        /*************************************************************************
        Reference FHT stub
        *************************************************************************/
        private static void RefFHT(ref double[] a, int n, int offs)
        {
            var buf = new double[n];
            var pin = Consts.pi2 / n;
            for (var i = 0; i < n; i++)
            {
                double v = 0;
                int j;
                for (j = 0; j <= n - 1; j++)
                {
                    var ij = i * j;
                    v += a[offs + j] * (Cos(pin * ij) + Sin(pin * ij));
                }
                buf[i] = v;
            }

            for (var i = 0; i <= n - 1; i++)
                a[offs + i] = buf[i];
        }

        #region Nested type: ftplan

        internal sealed class FT_Plan
        {
            public int[] plan;
            public double[] PreComputed;
            public double[] StackBuffer;
            public double[] TempBuffer;

            public FT_Plan()
            {
                plan = Array.Empty<int>();
                PreComputed = Array.Empty<double>();
                TempBuffer = Array.Empty<double>();
                StackBuffer = Array.Empty<double>();
            }
        }

        #endregion
    }

}