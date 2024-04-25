
MatrixZ.Test();

return;

const double fd = 2000;
const double dt = 1 / fd;

const double fp = 400;
const double fs = 300;
const double Rp = -1;
const double Rs = -40;

//const double df = 0.001;
//const int N = 10000;
//const double fd = N * df;
//const double dt = 1 / fd;
//var HH = Enumerable.Range(0, N)
//   .Select(i => i * df)
//   .Select(f => eq.GetTransmissionCoefficient(f, dt))
//   .ToAbs();

//var model = new PlotModel()
//   .SetBackground(OxyColors.White)
//   .Grid()
//   .Line(HH, df)
//   .SetMinY(0)
//   .ToPNG("test.png")
//   .Execute();

//return;

var Gp = Rp.From_dB();
var Gs = Rs.From_dB();


var filter = new ButterworthHighPass(dt, fs, fp, Gp, Gs);

var impulse_length = filter.GetImpulseResponse(Accuracy: Gs).Count();

var x = MathEnumerableSignal.Sin(dt, fs, (int)(50 / fp / dt));
var y = filter.ProcessIndividual(x);
var k = y.Power / x.Power;
var kdb = k.In_dB();
var x_a = x.PeakToPeakAmplitude;
var y_a = y.PeakToPeakAmplitude;
var k_a = y_a / x_a;

var Hf = filter.EnumTransmissionCoefficients(0, fd / 2, 1).ToArray();
var Hf_Abs = Hf.ToArray(v => v with { H = v.H.Power.In_dB_byPower() });
var Hf_Arg = Hf.ToArray(v => v with { H = v.H.Arg.ToDeg() });
var Ht = filter.EnumTransmissionCoefficientsTransient(dt, 20 / fp, 0, fd / 2, 1)
   .ToArray(v => v with { H = v.H.In_dB_byPower() });

const int count = 520;
//const double T = count * dt;
const double f_min = 295;
const double f_max = 305;

var file = new FileInfo($"test[{DateTime.Now:yyyy-MM-ddTHH-mm-ss}].png");
Signal.Random.SpectrumBand(dt, count, f_min, f_max)
   .Plot()
   .ToPNG(file)
   .Execute();


//await Task.Delay(3000);

//var process = Process.GetProcessesByName("Microsoft.Photos");

//if (process is null)
//{
//    Console.WriteLine("Complete!");
//    return;
//}

//var cancellation = new CancellationTokenSource();
//Console.CancelKeyPress += (_, e) =>
//{
//    Console.WriteLine("Ctrl+C pressed");
//    e.Cancel = true;
//    cancellation.Cancel();
//};

//await using var terminator = cancellation.Token.Register(p =>
//{
//    Console.WriteLine("Terminating window");
//    ((Process)p).CloseMainWindow();
//}, process);

//Console.WriteLine("Plot ready");

//try
//{
//    await process.WaitForExitAsync(cancellation.Token);

//    Console.WriteLine("Complete!");
//}
//catch (OperationCanceledException e) when(e.CancellationToken == cancellation.Token)
//{
//    Console.WriteLine("Operation terminated");
//}

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