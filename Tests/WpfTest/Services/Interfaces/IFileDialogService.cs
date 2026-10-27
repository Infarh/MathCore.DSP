#nullable enable
namespace WpfTest.Services.Interfaces;

/// <summary>Сервис выбора файлов в диалоговых окнах</summary>
public interface IFileDialogService
{
    /// <summary>Показать диалог открытия файла</summary>
    /// <param name="Title">Заголовок окна</param>
    /// <param name="Filter">Фильтр расширений</param>
    /// <param name="InitialDirectory">Начальный каталог</param>
    /// <returns>Выбранный путь или <c>null</c></returns>
    string? OpenFile(string Title, string Filter, string? InitialDirectory = null);

    /// <summary>Показать диалог сохранения файла</summary>
    /// <param name="Title">Заголовок окна</param>
    /// <param name="Filter">Фильтр расширений</param>
    /// <param name="DefaultFileName">Имя файла по умолчанию</param>
    /// <param name="InitialDirectory">Начальный каталог</param>
    /// <returns>Выбранный путь или <c>null</c></returns>
    string? SaveFile(string Title, string Filter, string DefaultFileName, string? InitialDirectory = null);
}
