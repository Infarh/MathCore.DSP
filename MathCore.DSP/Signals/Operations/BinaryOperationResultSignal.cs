namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал результата бинарной операции</summary>
public abstract class BinaryOperationResultSignal : OperationResultSignal
{
    /// <summary>Функция бинарной операции над отсчётами двух сигналов</summary>
    private readonly Func<double, double, double> _Function;

    /// <summary>Первый операнд-сигнал бинарной операции</summary>
    public DigitalSignal Signal1 { get; }

    /// <summary>Второй операнд-сигнал бинарной операции</summary>
    public DigitalSignal Signal2 { get; }

    /// <summary>Число отсчётов (вычисляется как число отсчётов первого сигнала)</summary>
    public override int SamplesCount => Signal1.SamplesCount;

    public override double this[int n]
    {
        get => _Function(Signal1[n], Signal2[n]);
        set => throw new NotSupportedException();
    }

    /// <summary>Инициализация нового цифрового сигнала бинарной операции</summary>
    /// <param name="Signal1">Первый сигнал-операнд операции</param>
    /// <param name="Signal2">Второй сигнал-операнд операции</param>
    /// <param name="Function">Функция, вычисляемая от каждой пары отсчётов сигналов</param>
    protected BinaryOperationResultSignal(
        DigitalSignal Signal1, 
        DigitalSignal Signal2, 
        Func<double, double, double> Function) 
        : base(Signal1.NotNull().dt, Math.Min(Signal1.t0, Signal2.NotNull().t0))
    {
        this.Signal1 = Signal1;
        this.Signal2 = Signal2;
        _Function = Function.NotNull();
    }

    public override IEnumerator<double> GetEnumerator() => Signal1.Zip(Signal2, _Function).GetEnumerator();
}