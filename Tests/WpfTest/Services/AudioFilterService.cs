#nullable enable
using System.IO;

using MathCore.DI;
using MathCore.DSP;
using MathCore.DSP.Filters;

using Microsoft.Extensions.DependencyInjection;

using NAudio.Wave;

namespace WpfTest.Services;

/// <summary>Аудиоданные, представленные межканальными отсчётами</summary>
/// <param name="Samples">Межканальные отсчёты в диапазоне -1 … +1</param>
/// <param name="SampleRate">Частота дискретизации, Гц</param>
/// <param name="Channels">Число каналов</param>
/// <param name="FileName">Имя исходного файла</param>
public sealed record AudioData(float[] Samples, int SampleRate, int Channels, string FileName)
{
    /// <summary>Число кадров (отсчётов на канал)</summary>
    public int FrameCount => Channels > 0 ? Samples.Length / Channels : 0;

    /// <summary>Длительность сигнала</summary>
    public TimeSpan Duration => TimeSpan.FromSeconds((double)FrameCount / SampleRate);

    /// <summary>Максимальное абсолютное значение отсчётов</summary>
    public double Peak
    {
        get
        {
            var peak = 0d;
            foreach (var sample in Samples)
            {
                var abs = Math.Abs(sample);
                if (abs > peak) peak = abs;
            }
            return peak;
        }
    }
}

/// <summary>Результат фильтрации аудиосигнала</summary>
/// <param name="Data">Отфильтрованные данные</param>
/// <param name="Gain">Коэффициент нормировки, применённый для предотвращения перегрузки</param>
public sealed record AudioFilterResult(AudioData Data, double Gain);

/// <summary>Загрузка, фильтрация, сохранение и воспроизведение аудиофайлов</summary>
[Service(ServiceLifetime.Singleton)]
public sealed class AudioFilterService : IDisposable
{
    private WaveOut? _Output;

    /// <summary>Загрузить аудиофайл (WAV, MP3, AIFF)</summary>
    /// <param name="Path">Путь к файлу</param>
    /// <returns>Загруженные аудиоданные</returns>
    public AudioData Load(string Path)
    {
        ArgumentNullException.ThrowIfNull(Path);
        if (!File.Exists(Path)) throw new FileNotFoundException("Аудиофайл не найден", Path);

        using var reader = new AudioFileReader(Path);
        var format = reader.WaveFormat;
        var samples = new List<float>((int)Math.Min(reader.Length / 4, 16_000_000));
        var buffer = new float[format.SampleRate * format.Channels];

        int read;
        while ((read = reader.Read(buffer.AsSpan())) > 0)
            samples.AddRange(buffer.AsSpan(0, read).ToArray());

        return new([.. samples], format.SampleRate, format.Channels, System.IO.Path.GetFileName(Path));
    }

    /// <summary>Отфильтровать аудиоданные с сохранением числа каналов</summary>
    /// <param name="Data">Исходные данные</param>
    /// <param name="Filter">Фильтр</param>
    /// <returns>Отфильтрованные данные</returns>
    public AudioFilterResult Filter(AudioData Data, IIR Filter)
    {
        ArgumentNullException.ThrowIfNull(Data);
        ArgumentNullException.ThrowIfNull(Filter);

        var a = Filter.A.ToArray();
        var b = Filter.B.ToArray();
        var samples = new float[Data.Samples.Length];
        var channels = Math.Max(1, Data.Channels);

        for (var channel = 0; channel < channels; channel++)
        {
            var state = new double[a.Length];
            for (var i = channel; i < Data.Samples.Length; i += channels)
                samples[i] = (float)state.FilterSample(a, b, Data.Samples[i]);
        }

        var peak = 0d;
        foreach (var sample in samples)
        {
            var abs = Math.Abs(sample);
            if (abs > peak) peak = abs;
        }

        var gain = peak > 0.99 ? 0.99 / peak : 1;
        if (gain < 1)
            for (var i = 0; i < samples.Length; i++) samples[i] = (float)(samples[i] * gain);

        return new(new(samples, Data.SampleRate, Data.Channels, Data.FileName), gain);
    }

    /// <summary>Сохранить данные в WAV-файл (16 бит PCM)</summary>
    /// <param name="Path">Путь к файлу</param>
    /// <param name="Data">Сохраняемые данные</param>
    public void Save(string Path, AudioData Data)
    {
        ArgumentNullException.ThrowIfNull(Path);
        ArgumentNullException.ThrowIfNull(Data);

        var bytes = ToPcm16(Data);
        using var writer = new WaveFileWriter(Path, new WaveFormat(Data.SampleRate, 16, Data.Channels));
        writer.Write(bytes, 0, bytes.Length);
    }

    /// <summary>Воспроизвести данные на устройстве вывода по умолчанию</summary>
    /// <param name="Data">Воспроизводимые данные</param>
    public void Play(AudioData Data)
    {
        ArgumentNullException.ThrowIfNull(Data);

        Stop();

        var bytes = ToPcm16(Data);
        var format = new WaveFormat(Data.SampleRate, 16, Data.Channels);
        var stream = new RawSourceWaveStream(bytes, 0, bytes.Length, format);
        var output = new WaveOut();
        output.PlaybackStopped += (_, _) =>
        {
            output.Dispose();
            stream.Dispose();
            if (ReferenceEquals(_Output, output)) _Output = null;
        };

        _Output = output;
        output.Init(stream);
        output.Play();
    }

    /// <summary>Остановить воспроизведение</summary>
    public void Stop()
    {
        var output = _Output;
        _Output = null;
        output?.Stop();
        output?.Dispose();
    }

    private static byte[] ToPcm16(AudioData Data)
    {
        var bytes = new byte[Data.Samples.Length * 2];
        for (var i = 0; i < Data.Samples.Length; i++)
        {
            var value = (short)(Math.Clamp(Data.Samples[i], -1f, 1f) * short.MaxValue);
            bytes[2 * i] = (byte)value;
            bytes[2 * i + 1] = (byte)(value >> 8);
        }
        return bytes;
    }

    /// <inheritdoc />
    public void Dispose() => Stop();
}
