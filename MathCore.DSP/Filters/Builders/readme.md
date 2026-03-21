# Fluent-интерфейс построения IIR-фильтров

Актуальный API построения IIR-фильтров расположен в `IirBuilder` и связанных типах.

## Точки входа

Есть три основных варианта старта:

- `Filter.IIR(double dt)` — старт с периодом дискретизации
- `Filter.IIRBySamplingFrequency(double fd)` — старт с частотой дискретизации
- `Filter.IIR().WithSampling(...)` / `Filter.IIR().WithSamplingFrequency(...)` — декларативный DSL-старт

## Поддерживаемые семейства и виды

Семейства:

- `Butterworth`
- `Chebyshev(...)`
- `Elliptic`
- `RC`
- `RLC`

Типы полос для IIR-семейств (`Butterworth/Chebyshev/Elliptic`):

- `LowPass`
- `HighPass`
- `BandPass`
- `BandStop`

## Базовые сценарии

### 1) Быстрое создание фильтра

```csharp
using MathCore.DSP.Filters;

var filter = Filter.IIR(1d / 10_000)
    .Butterworth
    .LowPass(1_000, 2_000)
    .Create();
```

### 2) Создание по `fd`

```csharp
using MathCore.DSP.Filters;

var filter = Filter.IIRBySamplingFrequency(10_000)
    .Elliptic
    .BandPass(500, 1_000, 2_000, 2_500)
    .WithGainsInDb(1, 40)
    .Create();
```

### 3) DSL-сценарий `WithPassband/WithStopband`

```csharp
using MathCore.DSP.Filters;

var filter = Filter.IIR()
    .WithSamplingFrequency(10_000)
    .Chebyshev(ChebyshevType.II)
    .BandStop()
    .WithPassband(500, 2_500)
    .WithStopband(1_000, 2_000)
    .BuildFilter();
```

## Работа со спецификацией

Интерфейс поддерживает раздельные этапы:

1. Сформировать спецификацию
2. Позже построить фильтр из спецификации

### 1) Получить спецификацию из builder

```csharp
using MathCore.DSP.Filters;
using MathCore.DSP.Filters.Builders;

IirFilterSpecification specification = Filter.IIR(1d / 10_000)
    .Butterworth
    .BandPass(500, 1_000, 2_000, 2_500)
    .BuildSpecification();
```

### 2) Собрать фильтр из спецификации

```csharp
using MathCore.DSP.Filters;
using MathCore.DSP.Filters.Builders;

Filter filter = specification.CreateFilter();
```

### 3) Безопасный Try-паттерн

```csharp
using MathCore.DSP.Filters;
using MathCore.DSP.Filters.Builders;

var success = Filter.IIR()
    .WithSamplingFrequency(10_000)
    .Butterworth
    .LowPass()
    .WithPassband(1_000)
    .TryGetSpecification(out var specification);

if (success)
{
    var filter = specification!.CreateFilter();
}
```

## Настройка коэффициентов передачи

Поддерживаются оба варианта:

- линейные коэффициенты: `WithGains(Gp, Gs)`
- значения в дБ: `WithGainsInDb(passRippleDb, stopAttenuationDb)`

Дополнительно:

- `WithPassbandRippleDb(...)`
- `WithStopbandAttenuationDb(...)`

## Иерархия спецификаций

Общая базовая абстракция:

- `IirFilterSpecification`

Вертикальная иерархия по семействам:

- `ButterworthSpecification`
- `ChebyshevSpecification`
- `EllipticSpecification`

Горизонтальные интерфейсы:

- `IIirSamplingSpecification`
- `IIirGainSpecification`
- `IIirFamilySpecification`
- `IIirPassTypeSpecification`
- `IIirLowHighFrequenciesSpecification`
- `IIirBandFrequenciesSpecification`
- `IChebyshevTypeSpecification`

Для простых веток также есть отдельные спецификации:

- `RcLowPassSpecification`
- `RcExponentialLowPassSpecification`
- `RcHighPassSpecification`
- `RlcBandPassSpecification`
- `RlcBandStopSpecification`

## RC/RLC примеры

```csharp
using MathCore.DSP.Filters;

var rc = Filter.IIR(1d / 1_000).RC.HighPass(50).Create();
var rlc = Filter.IIR(1d / 1_000).RLC.BandPass(120, 20).Create();
```

## Alias-методы

Для читаемости доступны алиасы:

- `BuildSpecification()` == `GetSpecification()`
- `BuildFilter()` == `Create()`
