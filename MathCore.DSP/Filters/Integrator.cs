namespace MathCore.DSP.Filters;

/// <summary>Интегратор</summary>
// Url: http://www.dsplib.ru/content/cic/cic.html
public class Integrator : IIR
{
    /// <summary>Инициализация нового интегратора</summary>
    public Integrator() : base([1d], [1d, -1d]) { }
}