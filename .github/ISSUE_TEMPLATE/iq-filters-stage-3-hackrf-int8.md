---
name: "IQ Filters / Этап 3 - HackRF int8 FIR/IIR"
about: "Доведение int8 ветки (HackRF One) до полного паритета"
title: "[IQ-Filters][Stage-3] HackRF int8 FIR/IIR parity"
labels: ["enhancement", "iq-filters", "stage-3", "hackrf"]
assignees: []
---

## Цель
Довести целочисленную int8 ветку HackRF до полного паритета FIR/IIR + stream API.

## Объём работ
- [ ] Довести FIR до паритета с уже реализованным IIR
- [ ] Унифицировать поведение ошибок и векторов состояния
- [ ] Проверить in-place/out-of-place, sync/async stream
- [ ] Зафиксировать квантование, округление и saturating clamp

## Критерии приёмки
- [ ] Полный набор unit-тестов для int8 покрывает FIR/IIR
- [ ] Совпадение с reференсом (complex/double + квантование)
- [ ] Тесты на формат interleaved и odd-length проходят

## Готовность к использованию
После этапа HackRF RX/TX можно использовать целиком в raw-представлении.

## Зависимости
- Stage-0
- Трекер: docs/IQFiltersTracking.md

## Checklist DoD
- [ ] Код
- [ ] Unit tests
- [ ] Документация
- [ ] Нет регрессий
