namespace MathCore.DSP.Signals.Operations
{
    public abstract class OperationResultSignal : DigitalSignal
    {
        protected OperationResultSignal(double dt) : base(dt) { }
    }
}