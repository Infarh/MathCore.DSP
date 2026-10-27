#nullable enable
using System.IO;

using Microsoft.Win32;

using WpfTest.Services.Interfaces;

namespace WpfTest.Services;

/// <summary>Сервис выбора файлов на основе стандартных диалогов WPF</summary>
public sealed class FileDialogService : IFileDialogService
{
    /// <inheritdoc />
    public string? OpenFile(string Title, string Filter, string? InitialDirectory = null)
    {
        var dialog = new OpenFileDialog
        {
            Title = Title,
            Filter = Filter,
            CheckFileExists = true,
        };

        if (Directory.Exists(InitialDirectory)) dialog.InitialDirectory = InitialDirectory;

        return dialog.ShowDialog() == true ? dialog.FileName : null;
    }

    /// <inheritdoc />
    public string? SaveFile(string Title, string Filter, string DefaultFileName, string? InitialDirectory = null)
    {
        var dialog = new SaveFileDialog
        {
            Title = Title,
            Filter = Filter,
            FileName = DefaultFileName,
            AddExtension = true,
            OverwritePrompt = true,
        };

        if (Directory.Exists(InitialDirectory)) dialog.InitialDirectory = InitialDirectory;

        return dialog.ShowDialog() == true ? dialog.FileName : null;
    }
}
