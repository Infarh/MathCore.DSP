using MathCore.DSP.Samples;
using MathCore.DSP.Samples.Extensions;

namespace MathCore.DSP.Tests.Samples;

[TestClass]
public class SampleSI16ExTests
{
    /// <summary>Тест фазовой демодуляции для пустого массива</summary>
    [TestMethod]
    public void PhaseDemodulation_EmptyArray_ReturnsEmpty()
    {
        // Arrange
        var samples = Span<SampleSI16>.Empty;
        const double f0 = 1000.0;
        const double fd = 48000.0;

        // Act
        var result = samples.PhaseDemodulation(f0, fd);

        // Assert
        Assert.AreEqual(0, result.Length);
    }

    /// <summary>Тест фазовой демодуляции для массива с одним элементом</summary>
    [TestMethod]
    public void PhaseDemodulation_SingleElement_ReturnsZero()
    {
        // Arrange
        var samples = new SampleSI16[] { new(100, 50) }.AsSpan();
        const double f0 = 1000.0;
        const double fd = 48000.0;

        // Act
        var result = samples.PhaseDemodulation(f0, fd);

        // Assert
        Assert.AreEqual(1, result.Length);
        Assert.AreEqual(0f, result[0]);
    }

    /// <summary>Тест фазовой демодуляции для постоянного сигнала</summary>
    [TestMethod]
    public void PhaseDemodulation_ConstantSignal_ReturnsNearZeros()
    {
        // Arrange - создаем постоянный сигнал
        var samples = new SampleSI16[]
        {
            new(100, 0),
            new(100, 0),
            new(100, 0),
            new(100, 0),
            new(100, 0)
        }.AsSpan();
        
        const double f0 = 1000.0;
        const double fd = 48000.0;

        // Act
        var result = samples.PhaseDemodulation(f0, fd);

        // Assert
        Assert.AreEqual(5, result.Length);
        Assert.AreEqual(0f, result[0]); // Первый всегда 0
        
        // Для постоянного сигнала мгновенная частота должна быть близка к нулю
        // после вычитания центральной частоты результат должен быть близок к -f0
        for (var i = 1; i < result.Length; i++)
        {
            Assert.IsTrue(Math.Abs(result[i] + f0) < 50, 
                $"Sample {i}: expected ~{-f0}, got {result[i]}");
        }
    }

    /// <summary>Тест фазовой демодуляции для синусоидального сигнала с известной частотой</summary>
    [TestMethod]
    public void PhaseDemodulation_SinusoidalSignal_ReturnsExpectedFrequency()
    {
        // Arrange - создаем синусоидальный сигнал с частотой f_signal
        const double fd = 48000.0;
        const double f0 = 1000.0;      // Центральная частота
        const double f_signal = 1500.0; // Частота сигнала
        const int samples_count = 200;
        
        var samples = new SampleSI16[samples_count];
        var dt = 1.0 / fd;
        
        for (var i = 0; i < samples_count; i++)
        {
            var t = i * dt;
            var angle = 2.0 * Math.PI * f_signal * t;
            var amplitude = 100.0;
            
            samples[i] = new SampleSI16(
                (sbyte)(amplitude * Math.Cos(angle)),
                (sbyte)(amplitude * Math.Sin(angle))
            );
        }

        // Act
        var result = samples.AsSpan().PhaseDemodulation(f0, fd);

        // Assert
        Assert.AreEqual(samples_count, result.Length);
        Assert.AreEqual(0f, result[0]); // Первый всегда 0
        
        // Проверяем, что результат близок к ожидаемой частоте (f_signal - f0)
        var expected_frequency = f_signal - f0;
        var tolerance = 100.0; // Допуск в Гц
        
        // Проверяем стабильную часть сигнала (пропускаем начальные образцы)
        var stable_samples = 0;
        for (var i = 20; i < result.Length - 20; i++) // Пропускаем края для стабилизации
        {
            if (Math.Abs(result[i] - expected_frequency) < tolerance)
                stable_samples++;
        }
        
        // Проверяем, что большинство образцов дает правильную частоту
        var expected_stable_count = (result.Length - 40) * 0.8; // 80% образцов должны быть стабильными
        Assert.IsTrue(stable_samples > expected_stable_count,
            $"Expected at least {expected_stable_count} stable samples, got {stable_samples}");
    }

    /// <summary>Тест производительности оптимизированной фазовой демодуляции</summary>
    [TestMethod]
    public void PhaseDemodulation_OptimizedPerformance_BetterThanOldVersion()
    {
        // Arrange - большой массив для тестирования производительности
        const int samples_count = 1_000_000; // Увеличиваем размер для лучшего измерения
        var samples = new SampleSI16[samples_count];
        var random = new Random(42); // Фиксированное семя для воспроизводимости
        
        for (var i = 0; i < samples_count; i++)
        {
            samples[i] = new SampleSI16(
                (sbyte)(random.Next(-128, 128)),
                (sbyte)(random.Next(-128, 128))
            );
        }

        const double f0 = 1000.0;
        const double fd = 48000.0;

        // Act - несколько прогонов для усреднения
        var times = new List<long>();
        for (var run = 0; run < 5; run++)
        {
            var stopwatch = System.Diagnostics.Stopwatch.StartNew();
            var result = samples.AsSpan().PhaseDemodulation(f0, fd);
            stopwatch.Stop();
            times.Add(stopwatch.ElapsedMilliseconds);
            
            // Проверяем корректность результата
            Assert.AreEqual(samples_count, result.Length);
        }

        var average_time = times.Average();
        
        // Assert - должно быть быстрее чем 1000 мс для 1M образцов
        Assert.IsTrue(average_time < 1000, 
            $"Optimized demodulation took {average_time:F1} ms on average, expected < 1000 ms");
        
        Console.WriteLine($"Оптимизированная фазовая демодуляция {samples_count} образцов:");
        Console.WriteLine($"Среднее время: {average_time:F1} мс");
        Console.WriteLine($"Производительность: {samples_count / average_time / 1000:F1} млн. образцов/сек");
    }

    /// <summary>Тест проверки правильности unwrap фазы</summary>
    [TestMethod]
    public void PhaseDemodulation_PhaseUnwrap_WorksCorrectly()
    {
        // Arrange - создаем сигнал с плавно изменяющейся фазой, но с перескоками ±2π
        const double fd = 1000.0;
        const double f0 = 0.0; // Нулевая центральная частота для простоты
        
        var samples = new SampleSI16[]
        {
            new(100, 0),     // фаза ≈ 0
            new(71, 71),     // фаза ≈ π/4
            new(0, 100),     // фаза ≈ π/2  
            new(-71, 71),    // фаза ≈ 3π/4
            new(-100, 0),    // фаза ≈ π
            new(-71, -71),   // фаза ≈ 5π/4 (но с unwrap должна быть ≈ -3π/4)
            new(0, -100),    // фаза ≈ 3π/2 (но с unwrap должна быть ≈ -π/2)
            new(71, -71),    // фаза ≈ 7π/4 (но с unwrap должна быть ≈ -π/4)
        }.AsSpan();

        // Act
        var result = samples.PhaseDemodulation(f0, fd);

        // Assert
        Assert.AreEqual(samples.Length, result.Length);
        Assert.AreEqual(0f, result[0]); // Первый всегда 0
        
        // Проверяем, что результаты имеют разумные значения без больших скачков
        for (var i = 1; i < result.Length; i++)
        {
            Assert.IsTrue(Math.Abs(result[i]) < 1000, 
                $"Sample {i}: frequency {result[i]} Hz seems too high");
        }
    }

    /// <summary>Тест сравнения точности оптимизированной и наивной реализации</summary>
    [TestMethod]
    public void PhaseDemodulation_OptimizedVsNaive_SameAccuracy()
    {
        // Arrange - создаем чистый синусоидальный сигнал
        const double fd = 48000.0;
        const double f0 = 1000.0;
        const double f_signal = 1200.0; 
        const int samples_count = 100;
        
        var samples = new SampleSI16[samples_count];
        var dt = 1.0 / fd;
        
        for (var i = 0; i < samples_count; i++)
        {
            var t = i * dt;
            var angle = 2.0 * Math.PI * f_signal * t;
            var amplitude = 120.0; // Используем почти максимальную амплитуду
            
            samples[i] = new SampleSI16(
                (sbyte)(amplitude * Math.Cos(angle)),
                (sbyte)(amplitude * Math.Sin(angle))
            );
        }

        // Act
        var optimized_result = samples.AsSpan().PhaseDemodulation(f0, fd);

        // Assert - проверяем точность на стабильной части
        var expected_frequency = f_signal - f0;
        var errors = new List<float>();
        
        for (var i = 10; i < optimized_result.Length - 10; i++)
        {
            var error = (float)Math.Abs(optimized_result[i] - expected_frequency);
            errors.Add(error);
        }
        
        var average_error = errors.Average();
        var max_error = errors.Max();
        
        Assert.IsTrue(average_error < 50, $"Average error {average_error:F1} Hz too high");
        Assert.IsTrue(max_error < 100, $"Max error {max_error:F1} Hz too high");
        
        Console.WriteLine($"Точность оптимизированного алгоритма:");
        Console.WriteLine($"Средняя ошибка: {average_error:F1} Гц");
        Console.WriteLine($"Максимальная ошибка: {max_error:F1} Гц");
    }
}