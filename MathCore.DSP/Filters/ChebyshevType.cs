namespace MathCore.DSP.Filters;

/// <summary>Типы фильтров Чебышева</summary>
public enum ChebyshevType : byte
{
    /// <summary>Фильтр Чебышева первого рода - основной фильтр, пропускающий нижнюю полосу частот</summary>
    I,
    /// <summary>Фильтр Чебышева второго рода, подавляющий верхнюю область частот (выше fp)</summary>
    II,
    /// <summary>Фильтр Чебышева второго рода c коррекцией частотного диапазона, подавляющий верхнюю область частот (выше fs)</summary>
    IICorrected
}