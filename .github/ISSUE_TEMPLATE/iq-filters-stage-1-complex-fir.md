---
name: "IQ Filters / Этап 1 - Complex FIR"
about: "Реализация комплексного FIR для IQ-потока"
title: "[IQ-Filters][Stage-1] Complex FIR"
labels: ["enhancement", "iq-filters", "stage-1"]
assignees: []
---

## Цель
Реализовать рабочий комплексный FIR для IQ с immediate value.

## Объём работ
- [ ] Реализовать `ComplexFIR` (состояние, Reset, Process)
- [ ] Добавить block API (`ReadOnlySpan<Complex> -> Span<Complex>`)
- [ ] Добавить/проверить `FrequencyResponse`
- [ ] Добавить потоковый API при необходимости

## Критерии приёмки
- [ ] Результаты совпадают с референсом (double по I/Q отдельно) в пределах eps
- [ ] Нет лишних аллокаций в hot-path Span API
- [ ] Тесты на состояние/сброс/границы входа проходят стабильно

## Готовность к использованию
После этапа можно сразу фильтровать complex IQ FIR в прод-пайплайне.

## Зависимости
- Stage-0
- Трекер: docs/IQFiltersTracking.md

## Checklist DoD
- [ ] Код
- [ ] Unit tests
- [ ] Документация с примером
- [ ] Нет регрессий
