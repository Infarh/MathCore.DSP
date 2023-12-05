using static System.Math;

// ReSharper disable UnusedType.Global
// ReSharper disable NotAccessedField.Local
// ReSharper disable InconsistentNaming

namespace MathCore.DSP.Fourier;

/// <summary>Элементарный фильтр Фурье</summary>
/// <remarks>Инициализация нового элементарного фильтра преобразования Фурье</remarks>
/// <param name="N">Размер выборки</param>
/// <param name="M">Размер спектра</param>
public class ElementaryFilter(int N, int M)
{
    /// <summary>Действительное значение</summary>
    private double _ReValue;
    /// <summary>Мнимое значение</summary>
    private double _ImValue;
#pragma warning restore IDE0052 // Удалить непрочитанные закрытые члены
    /// <summary>Дискрет фазы</summary>
    private readonly double _dPhi = M * Consts.pi2 / N;
    /// <summary>Индекс отсчёта</summary>
    private int _SamplesIndex;

    /// <summary>Значение фильтра</summary>
    public Complex Value => new(_ReValue / N, _ImValue / N);

    /// <summary>Инициализация фильтра (сброс состояния)</summary>
    public void Initialize()
    {
        _ReValue = 0;
        _ImValue = 0;
        _SamplesIndex = 0;
    }

    /// <summary>Обработка очередного значения</summary>
    /// <param name="value">Значение обрабатываемого отсчёта</param>
    public Complex Process(double value)
    {
        var arg = _dPhi * _SamplesIndex++;
        _ReValue += Cos(arg) * value;
        _ImValue += Sin(arg) * value;
        return Value;
    }
}