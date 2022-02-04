using System;
using System.Threading;
using System.Threading.Tasks;

namespace WpfTest.Infrastructure;

public sealed class ProcessTask<T> : IDisposable
{
    private readonly Action<IProgress<T>, CancellationToken> _Action;
    private readonly Action<Exception> _OnException;

    private readonly IProgress<T> _Progress;

    private readonly CancellationTokenSource _Cancellation = new();

    private Task _Task;

    public ProcessTask(
        Action<IProgress<T>, CancellationToken> Action, 
        Action<T> OnProgress = null,
        Action<Exception> OnException = null)
    {
        _Action = Action ?? throw new ArgumentNullException(nameof(Action));
        _OnException = OnException;
        _Progress = OnProgress is null ? null : new Progress<T>(OnProgress);
    }

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
            _OnException(e);
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