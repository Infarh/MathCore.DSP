using MathCore.DSP.Signals;

using static System.Math;

// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Filters;

/// <summary>Цифровой фильтр</summary>
public abstract class DigitalFilter : Filter
{
    /// <summary>Преобразование частоты аналогового  прототипа в частоту цифрового фильтра</summary>
    /// <param name="AnalogFrequency">Значение на оси частот цифрового фильтра</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Значение на оси частот аналогового прототипа</returns>
    public static double ToDigitalFrequency(double AnalogFrequency, double dt) => Tan(PI * AnalogFrequency * dt) / (PI * dt);

    /// <summary>Преобразование частоты аналогового  прототипа в частоту цифрового фильтра</summary>
    /// <param name="AnalogFrequencyW">Значение на оси частот цифрового фильтра</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Значение на оси частот аналогового прототипа</returns>
    public static double ToDigitalFrequencyW(double AnalogFrequencyW, double dt) => Tan(AnalogFrequencyW * dt / 2) / (dt / 2);

    /// <summary>Преобразование частоты цифрового фильтра в частоту аналогового прототипа</summary>
    /// <param name="DigitalFrequency">Значение на оси частот аналогового фильтра</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Значение на оси частот цифрового фильтра</returns>
    public static double ToAnalogFrequency(double DigitalFrequency, double dt) => Atan(PI * DigitalFrequency * dt) / (PI * dt);

    /// <summary>Преобразование частоты цифрового фильтра в частоту аналогового прототипа</summary>
    /// <param name="DigitalFrequencyW">Значение на оси частот аналогового фильтра</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Значение на оси частот цифрового фильтра</returns>
    public static double ToAnalogFrequencyW(double DigitalFrequencyW, double dt) => Atan(DigitalFrequencyW * dt / 2) / (dt / 2);

    /// <summary>Преобразование полюса из p-плоскости в z-плоскость</summary>
    /// <param name="p">Полюс p-плоскости</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Полюс в z-плоскости</returns>
    public static Complex ToZ(Complex p, double dt)
    {
        var w = 2 / dt;
        return (w + p) / (w - p);
    }

    /// <summary>Преобразование нулей/полюсов из p-плоскости в z-плоскость</summary>
    /// <param name="pPoles">Массив полюсов в p-плоскости</param>
    /// <param name="pZeros">Массив нулей в p-плоскости</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Массивы полюсов и нулей в z-плоскости</returns>
    /// <exception cref="ArgumentException">Если число нулей больше числа полюсов</exception>
    public static (Complex[] zPoles, Complex[] zZeros) ToZ(Complex[] pPoles, Complex[] pZeros, double dt)
    {
        if (pZeros.NotNull().Length > pPoles.NotNull().Length)
            throw new ArgumentException("Число нулей не должно превышать числа полюсов", nameof(pZeros));

        var poles_count = pPoles.Length;

        var z_zeros = new Complex[poles_count];
        var z_poles = new Complex[poles_count];

        for (var i = 0; i < pPoles.Length; i++) z_poles[i] = ToZ(pPoles[i], dt);
        for (var i = 0; i < pZeros.Length; i++) z_zeros[i] = ToZ(pZeros[i], dt);
        for (var i = pZeros.Length; i < z_zeros.Length; i++) z_zeros[i] = -1;

        return (z_poles, z_zeros);
    }

    /// <summary>Преобразование нулей/полюсов из p-плоскости в z-плоскость</summary>
    /// <param name="p">Перечисление нулей/полюсов p-плоскости</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Нули/полюса z-плоскости</returns>
    public static IEnumerable<Complex> ToZ(IEnumerable<Complex> p, double dt)
    {
        foreach (var z in p)
            yield return ToZ(z, dt);
    }

    /// <summary>Преобразование нулей/полюсов из p-плоскости в z-плоскость с масштабированием</summary>
    /// <param name="p">Перечисление нулей/полюсов p-плоскости</param>
    /// <param name="W0">Коэффициент масштабирования</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Нули/полюса z-плоскости с масштабированием</returns>
    public static IEnumerable<Complex> ToZ(IEnumerable<Complex> p, double W0, double dt)
    {
        foreach (var z in p)
            yield return ToZ(z * W0, dt);
    }

    /// <summary>Преобразование нулей/полюсов в массив в z-плоскости</summary>
    /// <param name="p">Нули/полюса p-плоскости</param>
    /// <param name="dt">Частота дискретизации</param>
    /// <param name="W0">Коэффициент масштабирования</param>
    /// <returns>Массив нулей/полюсов z-плоскости</returns>
    public static Complex[] ToZArray(IEnumerable<Complex> p, double dt, double W0 = 1) => ToZ(p, W0, dt).ToArray();

    /// <summary>Преобразование полюса из z-плоскости в p-плоскость</summary>
    /// <param name="z">Полюс z-плоскости</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Полюс в p-плоскости</returns>
    public static Complex ToP(Complex z, double dt) => 2 / dt * (z - 1) / (z + 1);

    /// <summary>Расчёт нормирующего множителя (приводящего системную-передаточную функцию к виду с максимумом в 1)</summary>
    /// <param name="poles">Набор полюсов</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Нормирующий множитель</returns>
    public static double GetNormalizeCoefficient(IEnumerable<Complex> poles, double dt)
    {
        if (poles is null) throw new ArgumentNullException(nameof(poles));

        var k = 2 / dt;
        var (re, im) = poles.Multiply(p => k - p);
        return (im / re).Abs() <= 1e-15
            ? 1 / re
            : throw new InvalidOperationException($"Вычисления привели к комплексному результату {new Complex(re, im)}");
    }

    /// <summary>Вектор состояния</summary>
    protected readonly double[] State;

    /// <summary>Порядок фильтра</summary>
    public virtual int Order => State.Length - 1;

    /// <summary>Инициализация нового цифрового фильтра</summary>
    /// <param name="Order">Порядок фильтра</param>
    /// <exception cref="ArgumentOutOfRangeException">Если порядок фильтра 0, или меньше</exception>
    protected DigitalFilter(int Order) => State = new double[Order > 0 ? Order : throw new ArgumentOutOfRangeException(nameof(Order), Order, "Порядок фильтра должен быть больше 0")];

    /// <summary>Обработать очередной отсчёт цифрового сигнала</summary>
    /// <param name="Sample">Обрабатываемый отсчёт цифрового сигнала</param>
    /// <param name="state">Вектор состояния фильтра</param>
    /// <returns>Значение сигнала на выходе фильтра после обработки отсчёта</returns>
    public abstract double Process(double Sample, double[] state);

    /// <summary>Обработать отсчёт цифрового сигнала</summary>
    /// <param name="Sample">Обрабатываемый отсчёт цифрового сигнала</param>
    /// <returns>Значение сигнала на выходе фильтра после обработки отсчёта</returns>
    public override double Process(double Sample) => Process(Sample, State);

    /// <summary>Обработать цифровой сигнал</summary>
    /// <param name="Signal">Цифровой сигнал</param>
    /// <param name="state">Вектор состояния фильтра</param>
    /// <returns>Обработанный цифровой сигнал</returns>

    public DigitalSignal Process(DigitalSignal Signal, double[] state)
    {
        if (state.NotNull().Length != Order + 1)
            throw new InvalidOperationException($"Длина вектора состояний {state.Length} не равна порядку фильтра {Order} + 1");

        return new SamplesDigitalSignal(Signal.NotNull().dt, Signal.Select(s => Process(s, state)));
    }

    /// <summary>Обработать цифровой сигнал независимо от состояния фильтра (вектор состояния создаётся на каждый вызов этого метода)</summary>
    /// <param name="Signal">Обрабатываемый цифровой сигнал</param>
    /// <returns>Обработанный цифровой сигнал</returns>
    public DigitalSignal ProcessIndividual(DigitalSignal Signal) => Process(Signal, new double[Order + 1]);

    /// <summary>Сбросить состояние фильтра</summary>
    public override void Reset() => Array.Clear(State, 0, State.Length);

    /// <summary>Получить коэффициент передачи фильтра на указанной частоте (КЧХ)</summary>
    /// <param name="f">Частота расчёта коэффициента передачи</param>
    /// <param name="dt">Период дискретизации</param>
    /// <returns>Комплексный коэффициент передачи фильтра</returns>
    public virtual Complex FrequencyResponse(double f, double dt) => FrequencyResponse(f * dt);
}