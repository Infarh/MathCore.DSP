using MathCore.DSP.Samples;
using MathCore.DSP.Samples.Extensions;

namespace MathCore.DSP.Examples;

/// <summary>Примеры использования фазовой модуляции и демодуляции</summary>
public static class PhaseModulationExamples
{
    /// <summary>Демонстрация базового круглого преобразования</summary>
    public static void BasicRoundTripExample()
    {
        Console.WriteLine("=== Базовый пример круглого преобразования ===");
        
        // Исходные данные для передачи (отклонения частоты в Гц)
        var original_data = new float[] { 0f, 200f, -300f, 150f, -100f, 0f };
        
        const double f0 = 2400e6;  // 2.4 ГГц центральная частота
        const double fd = 48000.0; // 48 кГц дискретизация
        
        Console.WriteLine("Исходные данные (отклонения частоты, Гц):");
        Console.WriteLine($"[{string.Join(", ", original_data.Select(x => x.ToString("F0")))}]");
        
        // Фазовая модуляция
        var modulated_samples = original_data.AsSpan().PhaseModulation(f0, fd, amplitude: 120f);
        Console.WriteLine($"\nМодулировано {modulated_samples.Length} I/Q отсчётов");
        
        // Фазовая демодуляция
        var demodulated_data = modulated_samples.AsSpan().PhaseDemodulation(f0, fd);
        Console.WriteLine($"Демодулировано {demodulated_data.Length} значений");
        
        // Сравнение результатов (пропускаем первые отсчёты из-за переходных процессов)
        Console.WriteLine("\nСравнение исходных и демодулированных данных:");
        Console.WriteLine("Индекс | Исходные | Демодул. | Ошибка");
        Console.WriteLine("-------|----------|----------|--------");
        
        for (var i = 2; i < Math.Min(original_data.Length, demodulated_data.Length); i++)
        {
            var error = Math.Abs(demodulated_data[i] - original_data[i]);
            Console.WriteLine($"{i,6} | {original_data[i],8:F1} | {demodulated_data[i],8:F1} | {error,6:F1}");
        }
    }
    
    /// <summary>Демонстрация непрерывной модуляции блоками</summary>
    public static void ContinuousBlockModulationExample()
    {
        Console.WriteLine("\n=== Пример непрерывной модуляции блоками ===");
        
        // Разбиваем данные на блоки
        var data_blocks = new[]
        {
            new float[] { 0f, 100f, 200f },
            new float[] { 300f, 200f, 100f },
            new float[] { 0f, -100f, -200f }
        };
        
        const double f0 = 915e6;    // 915 МГц ISM диапазон
        const double fd = 2e6;      // 2 МГц дискретизация
        
        var all_samples = new List<SampleSI16>();
        var phase = 0.0; // Начальная фаза
        
        for (var block_index = 0; block_index < data_blocks.Length; block_index++)
        {
            var data_block = data_blocks[block_index];
            Console.WriteLine($"\nБлок {block_index + 1}: [{string.Join(", ", data_block.Select(x => x.ToString("F0")))}] Гц");
            Console.WriteLine($"Начальная фаза: {phase:F3} рад");
            
            // Модуляция с сохранением фазы
            var (samples, final_phase) = data_block.AsSpan().PhaseModulation(f0, fd, phase);
            
            Console.WriteLine($"Финальная фаза: {final_phase:F3} рад");
            Console.WriteLine($"Сгенерировано {samples.Length} I/Q отсчётов");
            
            // Добавляем к общему потоку
            all_samples.AddRange(samples);
            
            // Продолжаем с сохранением фазы
            phase = final_phase;
        }
        
        Console.WriteLine($"\nВсего сгенерировано {all_samples.Count} I/Q отсчётов");
        
        // Демодуляция всего потока
        var demodulated = all_samples.ToArray().AsSpan().PhaseDemodulation(f0, fd);
        Console.WriteLine($"Демодулировано {demodulated.Length} значений");
        
        // Показываем восстановленные данные
        Console.WriteLine("\nВосстановленные данные (пропуская переходные процессы):");
        var flat_original = data_blocks.SelectMany(x => x).ToArray();
        
        for (var i = 3; i < Math.Min(flat_original.Length, demodulated.Length - 1); i++)
        {
            var error = Math.Abs(demodulated[i] - flat_original[i]);
            Console.WriteLine($"[{i}] Исходные: {flat_original[i]:F0} Гц, " +
                            $"Демодул.: {demodulated[i]:F1} Гц, " +
                            $"Ошибка: {error:F1} Гц");
        }
    }
    
    /// <summary>Тест производительности модуляции и демодуляции</summary>
    public static void PerformanceTest()
    {
        Console.WriteLine("\n=== Тест производительности ===");
        
        const int data_count = 1_000_000;
        const double f0 = 2.4e9;
        const double fd = 48000.0;
        
        // Генерируем тестовые данные
        var random = new Random(42);
        var test_data = new float[data_count];
        for (var i = 0; i < data_count; i++)
            test_data[i] = (float)(random.NextDouble() * 1000 - 500); // ±500 Гц
        
        Console.WriteLine($"Тестирование на {data_count:N0} образцах");
        
        // Тест модуляции
        var modulation_timer = System.Diagnostics.Stopwatch.StartNew();
        var modulated = test_data.AsSpan().PhaseModulation(f0, fd);
        modulation_timer.Stop();
        
        var modulation_rate = data_count / modulation_timer.Elapsed.TotalSeconds / 1e6;
        Console.WriteLine($"Модуляция: {modulation_timer.ElapsedMilliseconds} мс, " +
                         $"{modulation_rate:F1} млн. образцов/сек");
        
        // Тест демодуляции  
        var demodulation_timer = System.Diagnostics.Stopwatch.StartNew();
        var demodulated = modulated.AsSpan().PhaseDemodulation(f0, fd);
        demodulation_timer.Stop();
        
        var demodulation_rate = data_count / demodulation_timer.Elapsed.TotalSeconds / 1e6;
        Console.WriteLine($"Демодуляция: {demodulation_timer.ElapsedMilliseconds} мс, " +
                         $"{demodulation_rate:F1} млн. образцов/сек");
        
        // Проверка точности
        var errors = new List<float>();
        for (var i = 10; i < data_count - 10; i++) // Пропускаем края
        {
            var error = Math.Abs(demodulated[i] - test_data[i]);
            errors.Add(error);
        }
        
        var avg_error = errors.Average();
        var max_error = errors.Max();
        
        Console.WriteLine($"Средняя ошибка: {avg_error:F2} Гц");
        Console.WriteLine($"Максимальная ошибка: {max_error:F2} Гц");
        Console.WriteLine($"Точность: {100 * (1 - avg_error / 500):F1}%");
    }
    
    /// <summary>Запуск всех примеров</summary>
    public static void RunAllExamples()
    {
        Console.WriteLine("Примеры использования фазовой модуляции и демодуляции");
        Console.WriteLine("=====================================================");
        
        BasicRoundTripExample();
        ContinuousBlockModulationExample();
        PerformanceTest();
        
        Console.WriteLine("\n=== Примеры завершены ===");
    }
}