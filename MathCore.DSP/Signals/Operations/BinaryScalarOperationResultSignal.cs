namespace MathCore.DSP.Signals.Operations;

/// <summary>Сигнал-результат скалярной операции между сигналом и вещественным числом</summary>
public abstract class BinaryScalarOperationResultSignal : OperationResultSignal
{
    public delegate double SignalCombinator(double Sample, double value);

    /// <summary>Вещественное число, участвующее в операции с сигналом</summary>
    private readonly double _Value;

    /// <summary>Функция операции над отсчётами сигнала</summary>
    private readonly SignalCombinator _Function;

    /// <summary>Сигнал, участвующий в операции</summary>
    public DigitalSignal Signal { get; }

    /// <summary>Количество отсчётов</summary>
    public override int SamplesCount => Signal.SamplesCount;

    public override double this[int n]
    {
        get => _Function(Signal[n], _Value);
        set => throw new NotSupportedException();
    }

    /// <summary>Инициализация нового результирующего сигнала для сигнала и числа</summary>
    /// <param name="Signal">Исходный сигнал</param>
    /// <param name="Value">Вещественное число</param>
    /// <param name="Function">Функция, вызываемая для каждого отсчёта сигнала (первый аргумент) и вещественного числа (второй аргумент)</param>
    protected BinaryScalarOperationResultSignal(DigitalSignal Signal, double Value, SignalCombinator Function) : base(Signal.dt)
    {
        this.Signal = Signal ?? throw new ArgumentNullException(nameof(Signal));
        _Value = Value;
        _Function = Function ?? throw new ArgumentNullException(nameof(Function));
    }

    public override IEnumerator<double> GetEnumerator() => Signal.Select(s => _Function(s, _Value)).GetEnumerator();
}