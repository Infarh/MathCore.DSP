using OxyPlot;

namespace MathCore.DSP.Extensions;

public static class DigitalSignalEx
{
    public static PlotModel Plot(this DigitalSignal signal) =>
        new PlotModel()
           .Grid()
           .SetBackground(OxyColors.White)
           .Line(signal, OxyColors.Red, signal.dt, signal.t0);
}
