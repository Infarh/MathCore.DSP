using MathCore.DSP.Signals;
using MathCore.DSP.Signals.IO;

namespace MathCore.DSP;

class Program
{
    private static void RandomSignalsTest()
    {
        const double dt = 0.02;
        const double df = 0.04995004;
        const int ndf = 100;
        const double fd = 1 / dt;
        const double f0 = 0.1;
        const int count = 1001;

        var delta_f = new Interval(f0 - ndf * df, f0 + ndf * df);

        var random_signal = Signal.Random.SpectrumBand(dt, count, delta_f);

        var power = random_signal.Power;
    }

    static void Main()
    {
        //RandomSignalsTest();

        const string file_name = $"Sinus.signal16";
        const string path = $@"C:\123\tnd\Flight\TN2-Lightning1-operatori2iz2+zadatchiki-111.tnd.signals\{file_name} ";

        const string linear_signal_file = "linear.s16";

        const double dt = 0.1;
        const int samples_count = 1001;
        const double T = (samples_count - 1) * dt;
        const double A = T;
        var linear_signal = new SamplesDigitalSignal(dt, samples_count, t => A * t / T, -20);

        var linear_signal16 = new Signal16(linear_signal_file) { Max = A };
        linear_signal16.SetSignal(linear_signal);


        var signal = new Signal16(linear_signal_file);

        if (!signal.FileExists)
            throw new FileNotFoundException("Не найден файл данных сигнала", path);

        //var samples = signal.GetSamples(10, 15).ToArray();
        var samples = new List<SignalSample>(signal.SamplesCount);
        const double tmin = 10.123456789;
        const double t1 = 30, t2 = 40;
        foreach (var sample in signal.GetSamples(t1, t2))
        {
            samples.Add(sample);
        }
        samples.TrimExcess();

        var ssamples = linear_signal.Samples.Where(s => s.Time is >= t1 and <= t2).ToArray();

        var s0 = samples[0];
        var s1 = samples[^1];
        var ss0 = ssamples[0];
        var ss1 = ssamples[^1];

        if(s0.Time < tmin)
            Console.WriteLine(123);

        Console.WriteLine("  Total count: {0}", signal.Count());        // 23354
        Console.WriteLine("Samples count: {0}", samples.Count);        // 23354
        Console.WriteLine("          dt : {0}", signal.dt);             // 0,0009765625
        Console.WriteLine("          fd : {0}", 1 / signal.dt);         // 1024
        Console.WriteLine("          t0 : {0}", signal.t0);             // 0
        Console.WriteLine("          t1 : {0}", signal.t1);             // 22.806640625
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