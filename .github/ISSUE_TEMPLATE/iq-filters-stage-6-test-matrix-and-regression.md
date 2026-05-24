---
name: "IQ Filters / Этап 6 - Тестовая матрица и регрессия"
about: "Финализация тестовой системы для complex/int8/int16 IQ фильтров"
title: "[IQ-Filters][Stage-6] Полная тестовая матрица и регрессия"
labels: ["enhancement", "iq-filters", "stage-6", "tests"]
assignees: []
---

## Цель
Построить полную, стабильную и расширяемую тестовую матрицу для IQ-линейки.

## Объём работ
- [ ] Unit tests для Process/Reset/State
- [ ] Golden vectors для ключевых сценариев
- [ ] Property-based/инвариантные тесты
- [ ] Long-run/stability тесты
- [ ] Format tests (interleaved, endianness, odd-length, partial read)
- [ ] Cross-check complex vs double(I/Q), int vs quantized reference

## Критерии приёмки
- [ ] Покрытие критичных веток на целевом уровне
- [ ] Отсутствуют flaky тесты
- [ ] Есть профили fast/extended
- [ ] Все этапы 0-5 защищены регрессией

## Готовность к использованию
После этапа IQ-линейка готова к безопасной эволюции и рефакторингу.

## Зависимости
- Stages 0-5
- Трекер: docs/IQFiltersTracking.md

## Checklist DoD
- [ ] Код тестов
- [ ] CI профиль fast/extended
- [ ] Документация по тестовым профилям
- [ ] Нет регрессий
