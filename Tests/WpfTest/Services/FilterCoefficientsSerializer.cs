#nullable enable
using System.Globalization;
using System.IO;
using System.Text;
using System.Text.Encodings.Web;
using System.Text.Json;
using System.Text.Json.Serialization;

using MathCore.DI;

using Microsoft.Extensions.DependencyInjection;

using WpfTest.Models;

namespace WpfTest.Services;

/// <summary>Сохранение и загрузка коэффициентов БИХ-фильтра</summary>
/// <remarks>Поддерживается формат JSON и текстовый формат с секциями B и A</remarks>
[Service(ServiceLifetime.Singleton)]
public sealed class FilterCoefficientsSerializer
{
    private static readonly JsonSerializerOptions __JsonOptions = new()
    {
        WriteIndented = true,
        PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
        Encoder = JavaScriptEncoder.UnsafeRelaxedJsonEscaping,
        DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
    };

    /// <summary>Сериализовать коэффициенты в JSON</summary>
    /// <param name="Coefficients">Коэффициенты фильтра</param>
    /// <returns>Текст в формате JSON</returns>
    public string Serialize(FilterCoefficients Coefficients)
    {
        ArgumentNullException.ThrowIfNull(Coefficients);
        return JsonSerializer.Serialize(ToDto(Coefficients), __JsonOptions);
    }

    /// <summary>Загрузить коэффициенты из текста</summary>
    /// <param name="Text">Текст в формате JSON или в текстовом формате</param>
    /// <returns>Коэффициенты фильтра</returns>
    /// <exception cref="FormatException">Текст не содержит корректных коэффициентов</exception>
    public FilterCoefficients Deserialize(string Text)
    {
        if (string.IsNullOrWhiteSpace(Text))
            throw new FormatException("Файл коэффициентов пуст");

        return Text.TrimStart().StartsWith('{')
            ? FromDto(JsonSerializer.Deserialize<CoefficientsDto>(Text, __JsonOptions)
                      ?? throw new FormatException("Не удалось прочитать JSON коэффициентов"))
            : DeserializeText(Text);
    }

    /// <summary>Сохранить коэффициенты в файл</summary>
    /// <param name="Path">Путь к файлу</param>
    /// <param name="Coefficients">Коэффициенты фильтра</param>
    public void Save(string Path, FilterCoefficients Coefficients)
    {
        ArgumentNullException.ThrowIfNull(Path);
        ArgumentNullException.ThrowIfNull(Coefficients);
        File.WriteAllText(Path, Serialize(Coefficients), new UTF8Encoding(false));
    }

    /// <summary>Загрузить коэффициенты из файла</summary>
    /// <param name="Path">Путь к файлу</param>
    /// <returns>Коэффициенты фильтра</returns>
    public FilterCoefficients Load(string Path)
    {
        ArgumentNullException.ThrowIfNull(Path);
        if (!File.Exists(Path)) throw new FileNotFoundException("Файл коэффициентов не найден", Path);
        return Deserialize(File.ReadAllText(Path));
    }

    /// <summary>Текстовое представление коэффициентов для отображения</summary>
    /// <param name="Coefficients">Коэффициенты фильтра</param>
    /// <returns>Многострочный текст с полиномами B и A</returns>
    public string ToDisplayText(FilterCoefficients Coefficients)
    {
        ArgumentNullException.ThrowIfNull(Coefficients);
        var text = new StringBuilder();
        text.AppendLine($"B = [{string.Join(", ", Coefficients.B.Select(Format))}]");
        text.AppendLine($"A = [{string.Join(", ", Coefficients.A.Select(Format))}]");
        text.Append($"fd = {Coefficients.SamplingFrequency.ToString("F0", CultureInfo.InvariantCulture)} Гц");
        return text.ToString();

        static string Format(double Value) => Value.ToString("0.########", CultureInfo.InvariantCulture);
    }

    private static FilterCoefficients DeserializeText(string Text)
    {
        double[]? b = null;
        double[]? a = null;
        var sampling_frequency = 0d;

        foreach (var raw_line in Text.Split(['\r', '\n'], StringSplitOptions.RemoveEmptyEntries))
        {
            var line = raw_line.Trim();
            if (line.Length == 0) continue;

            var is_comment = line.StartsWith('#') || line.StartsWith("//");
            if (is_comment)
            {
                sampling_frequency = TryParseSamplingFrequency(line, sampling_frequency);
                continue;
            }

            var separator = line.IndexOfAny([':', '=']);
            if (separator > 0)
            {
                var key = line[..separator].Trim().ToLowerInvariant();
                var value = line[(separator + 1)..];
                switch (key)
                {
                    case "b":
                        b = ParseValues(value, "B");
                        continue;
                    case "a":
                        a = ParseValues(value, "A");
                        continue;
                    case "fd" or "fs" or "sampling" or "samplingfrequency" or "sample_rate" or "rate":
                        sampling_frequency = ParseOptionalFrequency(value, sampling_frequency);
                        continue;
                }
            }

            if (b is null) b = ParseValues(line, "B");
            else if (a is null) a = ParseValues(line, "A");
            else throw new FormatException("В файле обнаружено больше двух строк с коэффициентами");
        }

        if (b is null || a is null)
            throw new FormatException("В файле должны быть заданы коэффициенты B и A");

        return new(sampling_frequency > 0 ? sampling_frequency : 10_000, b, a, "Загружен из файла");
    }

    private static double TryParseSamplingFrequency(string Line, double Current)
    {
        if (Current > 0) return Current;
        var index = Line.IndexOfAny(['=', ':']);
        if (index < 0) return Current;
        var key = Line[..index].TrimStart('#', '/', ' ');
        if (!key.Contains("fd", StringComparison.OrdinalIgnoreCase)
            && !key.Contains("sampling", StringComparison.OrdinalIgnoreCase)
            && !key.Contains("rate", StringComparison.OrdinalIgnoreCase)
            && !key.Contains("частот", StringComparison.OrdinalIgnoreCase))
            return Current;
        var values = Line[(index + 1)..].Split([' ', '\t', ';'], StringSplitOptions.RemoveEmptyEntries);
        return values.Length > 0 && TryParseDouble(values[0], out var value) ? value : Current;
    }

    private static double ParseOptionalFrequency(string Text, double Current)
    {
        var value = Text.Split([' ', '\t', ';'], StringSplitOptions.RemoveEmptyEntries).FirstOrDefault();
        return value is not null && TryParseDouble(value, out var frequency) && frequency > 0 ? frequency : Current;
    }

    private static double[] ParseValues(string Text, string Name)
    {
        var tokens = Text.Split([' ', '\t', ';'], StringSplitOptions.RemoveEmptyEntries)
            .Select(v => v.Trim())
            .Where(v => v.Length > 0)
            .ToArray();

        // Одна строка с разделителем-запятой трактуется как CSV
        if (tokens.Length == 1 && tokens[0].Contains(','))
            tokens = tokens[0].Split(',', StringSplitOptions.RemoveEmptyEntries).Select(v => v.Trim()).ToArray();

        var values = tokens
            .Select(v => TryParseDouble(v, out var value)
                ? value
                : throw new FormatException($"Не удалось прочитать коэффициент «{v}» полинома {Name}"))
            .ToArray();

        return values.Length > 0 ? values : throw new FormatException($"Полином {Name} не содержит коэффициентов");
    }

    private static bool TryParseDouble(string Text, out double Value) =>
        double.TryParse(Text, NumberStyles.Float, CultureInfo.InvariantCulture, out Value)
        || double.TryParse(Text.Replace(',', '.'), NumberStyles.Float, CultureInfo.InvariantCulture, out Value);

    private static CoefficientsDto ToDto(FilterCoefficients Coefficients) => new()
    {
        Format = "MathCore.DSP.IIR",
        Version = 1,
        SamplingFrequency = Coefficients.SamplingFrequency,
        Description = Coefficients.Description,
        B = Coefficients.B,
        A = Coefficients.A,
    };

    private static FilterCoefficients FromDto(CoefficientsDto Dto)
    {
        if (Dto.A is null || Dto.A.Length < 2)
            throw new FormatException("В JSON отсутствует полином знаменателя A");

        if (Dto.B is null || Dto.B.Length == 0)
            throw new FormatException("В JSON отсутствует полином числителя B");

        if (Dto.B.Length > Dto.A.Length)
            throw new FormatException("Число коэффициентов B не должно превышать число коэффициентов A");

        return new(Dto.SamplingFrequency > 0 ? Dto.SamplingFrequency : 10_000, Dto.B, Dto.A, Dto.Description);
    }

    private sealed class CoefficientsDto
    {
        public string? Format { get; set; }
        public int Version { get; set; }
        public double SamplingFrequency { get; set; }
        public string? Description { get; set; }
        public double[]? B { get; set; }
        public double[]? A { get; set; }
    }
}
