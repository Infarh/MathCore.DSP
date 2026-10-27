#nullable enable
using System.Text;

using WpfTest.Models;
using WpfTest.Services;

namespace WpfTest.Infrastructure;

/// <summary>Самопроверка синтеза, анализа и сериализации фильтров</summary>
public static class SelfTest
{
    /// <summary>Выполнить самопроверку</summary>
    /// <param name="LogPath">Путь к файлу протокола</param>
    /// <returns>Код завершения: 0 — успех, 1 — ошибка</returns>
    public static int Run(string? LogPath = null)
    {
        var log = new StringBuilder();
        var errors = 0;
        var design = new FilterDesignService();
        var analysis = new FilterAnalysisService();
        var serializer = new FilterCoefficientsSerializer();

        void Write(string Message)
        {
            log.AppendLine(Message);
        }

        try
        {
            foreach (var kind_info in FilterKindInfo.All)
            {
                foreach (var family in FilterFamilyInfo.All.Where(f => f.Supports(kind_info.Kind)))
                {
                    var settings = CreateSettings(kind_info.Kind, family.Family);
                    FilterCoefficients? custom = null;

                    if (family.Family == FilterFamily.Custom)
                    {
                        var base_settings = CreateSettings(kind_info.Kind, FilterFamily.Butterworth);
                        var base_result = design.Create(base_settings);
                        custom = base_result.GetCoefficients();
                    }

                    try
                    {
                        var result = design.Create(settings, custom);
                        var result_analysis = analysis.Calculate(result.Filter, settings.SamplingFrequency);

                        Check(result.Order > 0, $"{kind_info.ShortTitle}/{family.Family}: порядок {result.Order}");
                        Check(result_analysis.Frequencies.Length == 2048, "число точек АЧХ");
                        Check(result_analysis.MagnitudeDb.All(double.IsFinite), "конечность значений АЧХ");
                        Check(result_analysis.PhaseDegrees.All(v => !double.IsInfinity(v)), "конечность значений ФЧХ");
                        Check(result_analysis.ImpulseResponse.Length > 0, "импульсная характеристика");
                        Check(result_analysis.Poles.Length == result.Filter.A.Count - 1, "число полюсов");
                        Check(result_analysis.Zeros.Length <= result.Order, "число нулей");
                        Check(result_analysis.Poles.All(p => p.Magnitude < 1 + 1e-3), "устойчивость фильтра");

                        Write($"OK  {kind_info.ShortTitle,-4} {family.Family,-22} порядок={result.Order,-3} B={result.Filter.B.Count} A={result.Filter.A.Count} max|p|={result_analysis.Poles.Max(p => p.Magnitude):F9}");
                    }
                    catch (Exception error)
                    {
                        errors++;
                        Write($"ERR {kind_info.ShortTitle,-4} {family.Family,-22} {error.Message}");
                    }
                }
            }

            var coefficients = design.Create(CreateSettings(FilterKind.BandPass, FilterFamily.Butterworth)).GetCoefficients();
            var json = serializer.Serialize(coefficients);
            var json_copy = serializer.Deserialize(json);
            Check(AreEqual(coefficients.B, json_copy.B) && AreEqual(coefficients.A, json_copy.A), "JSON-сериализация");

            var text = $"B: {string.Join(' ', coefficients.B)}\r\nA: {string.Join(' ', coefficients.A)}";
            var text_copy = serializer.Deserialize(text);
            Check(AreEqual(coefficients.A, text_copy.A), "текстовый формат");

            var temp_path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), $"mathcore_dsp_coefficients_{Guid.NewGuid():N}.json");
            try
            {
                serializer.Save(temp_path, coefficients);
                var file_copy = serializer.Load(temp_path);
                Check(AreEqual(coefficients.B, file_copy.B), "сохранение/загрузка файла");
            }
            finally
            {
                if (System.IO.File.Exists(temp_path)) System.IO.File.Delete(temp_path);
            }

            Write("OK  Сериализация коэффициентов (JSON, текст, файл)");

            var audio = new AudioFilterService();
            const int rate = 10_000;
            const int samples_count = rate;
            var samples = new float[samples_count];
            for (var i = 0; i < samples_count; i++)
                samples[i] = (float)(0.4 * Math.Sin(2 * Math.PI * 440 * i / rate) + 0.4 * Math.Sin(2 * Math.PI * 4_000 * i / rate));

            var audio_path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), $"mathcore_dsp_audio_{Guid.NewGuid():N}.wav");
            try
            {
                audio.Save(audio_path, new(samples, rate, 1, "selftest.wav"));
                var loaded = audio.Load(audio_path);
                Check(loaded.SampleRate == rate, "частота дискретизации WAV");
                Check(Math.Abs(loaded.FrameCount - samples_count) <= 1, "длина WAV");
                Check(loaded.Channels == 1, "число каналов WAV");

                var low_pass = design.Create(CreateSettings(FilterKind.LowPass, FilterFamily.Butterworth)).Filter;
                var filtered = audio.Filter(loaded, low_pass).Data;
                var amplitude_440 = GetAmplitude(filtered.Samples, rate, 440);
                var amplitude_4000 = GetAmplitude(filtered.Samples, rate, 4_000);

                Check(amplitude_440 > 0.25, $"прохождение 440 Гц (амплитуда {amplitude_440:F3})");
                Check(amplitude_4000 < 0.02, $"подавление 4000 Гц (амплитуда {amplitude_4000:F4})");

                Write($"OK  Аудио: WAV сохранение/загрузка и фильтрация (440 Гц={amplitude_440:F3}, 4 кГц={amplitude_4000:F4})");
            }
            finally
            {
                if (System.IO.File.Exists(audio_path)) System.IO.File.Delete(audio_path);
            }
        }
        catch (Exception error)
        {
            errors++;
            Write($"ERR {error}");
        }

        var summary = errors == 0 ? "САМОПРОВЕРКА ПРОЙДЕНА" : $"ОБНАРУЖЕНО ОШИБОК: {errors}";
        Write(summary);

        if (LogPath is not null)
        {
            try
            {
                System.IO.File.WriteAllText(LogPath, log.ToString(), new UTF8Encoding(false));
            }
            catch
            {
                // Протокол не критичен для результата проверки
            }
        }

        return errors == 0 ? 0 : 1;

        static void Check(bool Condition, string Message)
        {
            if (!Condition) throw new InvalidOperationException($"Проверка не пройдена: {Message}");
        }

        static bool AreEqual(double[] A, double[] B) => A.Length == B.Length && A.Zip(B).All(p => Math.Abs(p.First - p.Second) < 1e-12);

        static double GetAmplitude(float[] Samples, int Rate, double Frequency)
        {
            var (re, im) = (0d, 0d);
            for (var i = 0; i < Samples.Length; i++)
            {
                var phase = 2 * Math.PI * Frequency * i / Rate;
                re += Samples[i] * Math.Cos(phase);
                im -= Samples[i] * Math.Sin(phase);
            }
            return 2 * Math.Sqrt(re * re + im * im) / Samples.Length;
        }
    }

    private static FilterSettings CreateSettings(FilterKind Kind, FilterFamily Family)
    {
        var settings = new FilterSettings { Kind = Kind, Family = Family, SamplingFrequency = 10_000 };

        switch (Kind)
        {
            case FilterKind.LowPass:
                settings.PassLowFrequency = 1_000;
                settings.StopLowFrequency = 2_000;
                break;

            case FilterKind.HighPass:
                settings.StopLowFrequency = 500;
                settings.PassLowFrequency = 1_000;
                break;

            case FilterKind.BandPass:
                settings.StopLowFrequency = 500;
                settings.PassLowFrequency = 1_000;
                settings.PassHighFrequency = 2_000;
                settings.StopHighFrequency = 2_500;
                break;

            case FilterKind.BandStop:
                settings.PassLowFrequency = 500;
                settings.StopLowFrequency = 1_000;
                settings.StopHighFrequency = 2_000;
                settings.PassHighFrequency = 2_500;
                break;
        }

        return settings;
    }
}
