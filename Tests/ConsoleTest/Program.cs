using System;
using System.Linq.Expressions;

namespace MathCore.DSP
{
    class Program
    {
        static void Main()
        {
            static void Print(Expression<Func<double, double>> f) => Console.WriteLine(f);

            Print(x => Math.Sin(x) / x);
        }
    }
}
