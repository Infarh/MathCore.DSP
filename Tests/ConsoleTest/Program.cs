using System;
using System.Linq;

namespace MathCore.DSP
{
    class Program
    {
        static void Main()
        {
            var X = Enumerable.Range(1, 10).ToArray();

            Array.Copy(X, 0, X, 1, X.Length - 1);


            Console.ReadLine();
        }
    }
}
