using MathCore.DSP.Extensions;
using MathCore.DSP.Signals;
using MathCore.DSP.Signals.IO;

using OxyPlot;

namespace MathCore.DSP;

class Program
{
    static void Main()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const int count = 520;
        const double T = count * dt;
        const double f_min = 295;
        const double f_max = 305;

        Signal.Random.SpectrumBand(dt, count, f_min, f_max)
           .Plot()
           .ToPNG($"test[{DateTime.Now:yyyy-MM-ddTHH-mm-ss}].png")
           .Execute();
    }
}

//internal static class Ext
//{
//    public static double Median(this ICollection<double> items)
//    {
//        var items_count = items.Count;
//        switch (items_count)
//        {
//            case < 2: return double.NaN;
//            case 2: return items.Average();
//        }

//        var mediane = items.Average();
//        //var mm = new List<double>() { mediane };
//        int delta;
//        do
//        {
//            var count = 0;
//            foreach (var item in items)
//                if (item > mediane)
//                    count++;

//            delta = 2 * count - items_count;
//            var accuracy = delta / (double)items_count;
//            mediane += mediane * accuracy;
//            //mm.Add(mediane);
//        }
//        while (Math.Abs(delta) > 1);

//        return mediane;
//    }
//}