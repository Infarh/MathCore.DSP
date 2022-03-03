using System.ComponentModel.DataAnnotations;
using System.Security.Cryptography;

using MathCore.DSP.Infrastructure;
using MathCore.DSP.Signals;
using MathCore.DSP.Signals.IO;
using MathCore.Statistic.RandomNumbers;

namespace MathCore.DSP;

class Program
{
    private static void RandomSignalsTest()
    {
        const double fd = 1000;
        const double dt = 1 / fd;
        const int count = 1000;
        const double df = fd / count;
        const double f_min = 300;
        const double f_max = 400;

        var rnd = new Random();

        var random_signal = Signal.Random.SpectrumBand(dt, count, f_min, f_max, rnd);

        var power = random_signal.Power;
    }

    static void Main()
    {
        RandomSignalsTest();

        const int samples_count = 10000;
        var xx = Enumerable.Range(0, samples_count + 1);
        const double dt = 0.1;



        static double Signal(double t) => Math.Sin(2 * Math.PI * t);

        var rnd = new Random(5);
        var ss = xx.Select(t => /*10 * Signal(t * dt) +*/ (rnd.NextDouble() - 0.5) * 2);

        const string signal_file = "sin.signal16";

        //if (!File.Exists(signal_file))
        {
            var signal16 = new Signal16(signal_file) { t0 = 5, dt = dt, Min = -10, Max = 10 };
            signal16.SetSamples(ss);
        }

        var y = new Signal16(signal_file);

        var samples = y.GetSamples().ToArray();

        var pack = SamplesPack.Create(samples.AsDouble(), 10);

        var X = rnd.NextUniform(10000000, 17, Math.PI);
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