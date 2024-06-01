namespace WpfTest.Infrastructure;

public sealed class ProcessTask<T>(
    Action<IProgress<T>, CancellationToken> Action,
    Action<T> OnProgress = null,
    Action<Exception> OnException = null)
    : IDisposable
{
    private readonly Action<IProgress<T>, CancellationToken> _Action = Action ?? throw new ArgumentNullException(nameof(Action));

    private readonly IProgress<T> _Progress = OnProgress is null ? null : new Progress<T>(OnProgress);

    private readonly CancellationTokenSource _Cancellation = new();

    private Task _Task;

    public async void Start()
    {
        if (_Task is not null) return;
        lock (_Cancellation)
        {
            if (_Task is not null) return;
            _Task = Task.Run(TaskAction);
        }

        try
        {
            await _Task;
        }
        catch (OperationCanceledException) { }
        catch (Exception e)
        {
            OnException(e);
        }
        _Task = null;
    }

    public void Cancel() => _Cancellation.Cancel();

    private void TaskAction() => _Action(_Progress, _Cancellation.Token);

    public void Dispose()
    {
        _Cancellation.Cancel();
        _Cancellation.Dispose();
    }
}