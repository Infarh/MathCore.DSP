namespace MathCore.DSP.Filters.Builders;

/// <summary>Построитель полосно-задерживающего фильтра</summary>
/// <param name="dt">Период дискретизации</param>
public readonly record struct BandStopBuilder(double dt);
