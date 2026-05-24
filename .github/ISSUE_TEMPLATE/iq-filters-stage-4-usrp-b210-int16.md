---
name: "IQ Filters / Этап 4 - USRP B210 int16 LE"
about: "Реализация FIR/IIR для interleaved int16 LE IQ потока"
title: "[IQ-Filters][Stage-4] USRP B210 int16 LE FIR/IIR"
labels: ["enhancement", "iq-filters", "stage-4", "usrp-b210"]
assignees: []
---

## Цель
Реализовать полноценную int16 LE ветку для USRP B210 (FIR/IIR + stream API).

## Объём работ
- [ ] Добавить обработку interleaved int16 LE без лишних копирований
- [ ] Реализовать FIR/IIR block API (in-place/out-of-place)
- [ ] Добавить stream sync/async методы
- [ ] Проверить endian/sign/range поведение

## Критерии приёмки
- [ ] Паритет API с веткой HackRF int8
- [ ] Совпадение с референсом на эталонных наборах
- [ ] Тесты на endian и граничные значения проходят

## Готовность к использованию
После этапа можно подключать USRP B210 pipeline без адаптеров вне библиотеки.

## Зависимости
- Stage-0
- Stage-2/3 (желательно)
- Трекер: docs/IQFiltersTracking.md

## Checklist DoD
- [ ] Код
- [ ] Unit tests
- [ ] Документация
- [ ] Нет регрессий
