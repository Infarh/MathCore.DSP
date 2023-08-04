namespace MathCore.DSP.Signals;

public class SamplesPack
{
    public readonly struct Value(double Mu, double Sgm, double Min, double Max)
    {
        public double Mu { get; } = Mu;
        public double Sgm { get; } = Sgm;
        public double Min { get; } = Min;
        public double Max { get; } = Max;

        public double Delta => Max - Min;

        public double MinSgm1 => Mu - Sgm;
        public double MinSgm2 => Mu - Sgm * 2;
        public double MinSgm3 => Mu - Sgm * 3;
        public double MaxSgm1 => Mu + Sgm;
        public double MaxSgm2 => Mu + Sgm * 2;
        public double MaxSgm3 => Mu + Sgm * 3;

        public Interval MinMax => new(Min, Max);
        public double MinMaxMiddle => (Min + Max) / 2;
        public Interval Sgm1 => new(MinSgm1, MaxSgm1);
        public Interval Sgm2 => new(MinSgm2, MaxSgm2);
        public Interval Sgm3 => new(MinSgm3, MaxSgm3);

        public override string ToString() => $"{Mu}(±{3 * Sgm})[{Min}..{Max}]";
    }

    public static SamplesPack Create(IEnumerable<double> Samples, int PackSize)
    {
        var items = new List<Value>();

        var min = double.PositiveInfinity;
        var max = double.NegativeInfinity;
        var sum = 0d;
        var sum2 = 0d;
        var count = 0;
        foreach (var sample in Samples)
            if (count == PackSize)
            {
                var mu = sum / PackSize;
                var sgm = Math.Sqrt(sum2 / PackSize - mu * mu);
                items.Add(new(mu, sgm, min, max));

                min = sample;
                max = sample;
                sum = sample;
                sum2 = sample * sample;
                count = 1;
            }
            else
            {
                sum += sample;
                sum2 += sample * sample;
                min = Math.Min(min, sample);
                max = Math.Max(max, sample);
                count++;
            }

        if (count > 0)
        {
            var mu = sum / PackSize;
            var sgm = Math.Sqrt(sum2 - mu * mu);
            items.Add(new(mu, sgm, min, max));
        }

        return new(items.ToArray());
    }

    public IReadOnlyList<Value> Values { get; }

    public IEnumerable<double> Mu => Values.Select(e => e.Mu);

    public IEnumerable<double> Min => Values.Select(e => e.Min);
    public IEnumerable<double> Max => Values.Select(e => e.Max);
    public IEnumerable<double> Sgm => Values.Select(e => e.Sgm);

    public IEnumerable<double> MinSgm1 => Values.Select(e => e.MinSgm1);
    public IEnumerable<double> MinSgm2 => Values.Select(e => e.MinSgm2);
    public IEnumerable<double> MinSgm3 => Values.Select(e => e.MinSgm3);
    public IEnumerable<double> MaxSgm1 => Values.Select(e => e.MaxSgm1);
    public IEnumerable<double> MaxSgm2 => Values.Select(e => e.MaxSgm2);
    public IEnumerable<double> MaxSgm3 => Values.Select(e => e.MaxSgm3);

    public IEnumerable<Interval> MinMax => Values.Select(v => v.MinMax);
    public IEnumerable<double> MinMaxMiddle => Values.Select(v => v.MinMaxMiddle);
    public IEnumerable<Interval> Sgm1 => Values.Select(v => v.Sgm1);
    public IEnumerable<Interval> Sgm2 => Values.Select(v => v.Sgm2);
    public IEnumerable<Interval> Sgm3 => Values.Select(v => v.Sgm3);

    private SamplesPack(IReadOnlyList<Value> Values) => this.Values = Values;
}