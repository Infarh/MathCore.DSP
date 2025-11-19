# MathCore.DSP

Библиотека алгоритмов цифровой обработки сигналов (DSP) для .NET.
Предоставляет набор классов и функций для спектрального анализа, проектирования и применения цифровых фильтров (IIR/FIR), преобразований Фурье, оконных функций и вспомогательных операций над сигналами.

## Возможности

- Быстрое преобразование Фурье (FFT) для действительных и комплексных сигналов
- Реализация алгоритмов Bluestein и Rader для произвольной/простой длины последовательности
- Проектирование фильтров нижних/верхних/полосовых частот (Butterworth, Chebyshev, Elliptic, RC и др.)
- FIR‑фильтры и базовые операции свёртки
- Набор оконных функций: Rectangular, Hann, Hamming, Blackman, Blackman–Harris, Blackman–Nuttall, Nuttall, Flat Top, Bartlett–Hann и др.
- Вычисление АЧХ фильтров, частотный отклик, частотные характеристики
- Утилиты для анализа спектров, работы с комплексными числами, ресэмплинг

## Установка

Добавьте пакет через менеджер пакетов NuGet:
```
Install-Package MathCore.DSP
```
или через .NET CLI:
```
dotnet add package MathCore.DSP
```

## Быстрый старт

### FFT
```csharp
using MathCore.DSP.Fourier;

var signal = Enumerable.Range(0, 1024)
    .Select(i => Math.Sin(2 * Math.PI * 50 * i / 1024))
    .ToArray();

var spectrum = signal.FastFourierTransform(); // Комплексный спектр
```

### Применение окна
```csharp
using MathCore.DSP.WindowFunctions;

int N = 1024;
var windowed = new double[N];
for (var n = 0; n < N; n++)
    windowed[n] = HannWindow.Value(n, N) * signal[n];
```

### Проектирование фильтра Баттерворта НЧ
```csharp
using MathCore.DSP.Filters;

double fd = 10_000;            // Частота дискретизации
double dt = 1 / fd;            // Период дискретизации

double fp = 800;               // Граничная частота полосы пропускания (Гц)
double fs = 2000;              // Граничная частота полосы подавления (Гц)

var filter = new ButterworthLowPass(dt, fp, fs); // Gp=-1дБ, Gs=-40дБ по умолчанию

// Фильтрация потока отсчётов
foreach (var x in signal)
{
    var y = filter.Process(x); // y - отфильтрованный отсчёт
}
```

### Использование строителей фильтров
```csharp
using MathCore.DSP.Filters.Builders;

double fd = 48_000;
var lp = new LowPassBuilder(1 / fd)
    .Butterworth(fs: 6000, fp: 3000, Gp: 0.891250938, Gs: 0.01) // Параметры
    .Create();
```

### FIR фильтр
```csharp
using MathCore.DSP.Filters;

// Пример простейшего усредняющего FIR
var h = Enumerable.Repeat(1.0 / 8, 8).ToArray();
var fir = new FIR(h);
var y_array = signal.Select(fir.Process).ToArray();
```

## Основные классы

- `FFT` / `fft` – реализации алгоритмов БПФ
- `ButterworthLowPass`, `ChebyshevLowPass`, `EllipticLowPass`, `RCLowPass` – IIR‑фильтры
- `FIR` – фильтр с конечной импульсной характеристикой
- `LowPassBuilder`, `HighPassBuilder`, `BandPassBuilder` – строители для удобного создания фильтров (если доступны в сборке)
- Оконные функции в пространстве имён `MathCore.DSP.WindowFunctions.*`

## Производительность

- Реализованы оптимизации для степеней 2, а также Bluestein/Rader для произвольных длин
- Возможна асинхронная рекурсивная версия FFT (`Recursive_FFTAsync`) для больших массивов

## Лицензия

MIT. Свободно для использования в коммерческих и некоммерческих проектах при сохранении копирайта.

## Обратная связь

- Репозиторий: https://github.com/Infarh/MathCore.DSP
- Issues / предложения: через GitHub Issues

## Примечания

Некоторые функции могут требовать пакет `MathCore` (уже указан в зависимостях). Для .NET Standard 2.0 могут быть недоступны отдельные современные языковые конструкции.

Если вам нужны дополнительные типы окон или проектирование FIR по частотным маскам – создайте запрос.
