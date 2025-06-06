namespace MathCore.DSP.Filters.Builders;

/// <summary>Построитель полосно-пропускающего фильтра</summary>
/// <param name="dt">Период дискретизации</param>
public readonly record struct BandPassBuilder(double dt);
