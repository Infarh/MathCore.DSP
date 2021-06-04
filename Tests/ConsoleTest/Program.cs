using System;
using System.Linq.Expressions;

using static System.Math;

namespace MathCore.DSP
{
    class Program
    {
        static void Main()
        {
            static void Print(Expression<Func<double, double>> f) => Console.WriteLine(f);

            Print(x => Sin(x) / x);
        }
    }
}
