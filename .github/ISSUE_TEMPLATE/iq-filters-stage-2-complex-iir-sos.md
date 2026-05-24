---
name: "IQ Filters / Этап 2 - Complex IIR + SOS"
about: "Реализация комплексного IIR с поддержкой SOS"
title: "[IQ-Filters][Stage-2] Complex IIR + SOS"
labels: ["enhancement", "iq-filters", "stage-2"]
assignees: []
---

## Цель
Реализовать complex IIR для IQ-потока, включая устойчивый путь через SOS.

## Объём работ
- [ ] Реализовать `ComplexIIR` (direct form)
- [ ] Реализовать SOS-путь для высоких порядков
- [ ] Повторить методику валидации секций
- [ ] Обеспечить parity с существующим `IIR`

## Критерии приёмки
- [ ] Стабильность на длинных прогонах подтверждена тестами
- [ ] Сходимость с референсом (double I/Q) в пределах заданного допуска
- [ ] Тесты на частотный отклик и корректность состояния проходят

## Готовность к использованию
После этапа доступны production-ready complex IIR фильтры.

## Зависимости
- Stage-0
- Stage-1 (желательно)
- Трекер: docs/IQFiltersTracking.md

## Checklist DoD
- [ ] Код
- [ ] Unit tests + long-run tests
- [ ] Документация
- [ ] Нет регрессий
