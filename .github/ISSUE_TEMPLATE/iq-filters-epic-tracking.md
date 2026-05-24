---
name: "IQ Filters / Epic tracking"
about: "Эпик-трекер всей программы работ по IQ-фильтрам"
title: "[IQ-Filters][Epic] Комплексные и целочисленные IQ-фильтры"
labels: ["tracking", "iq-filters", "epic"]
assignees: []
---

## Эпик
Единая программа внедрения IQ-фильтров для complex, HackRF int8 и USRP B210 int16 LE.

## Источник плана
- docs/IQFiltersTracking.md

## Этапы
- [ ] Stage-0: Контракты и базовые адаптеры
- [ ] Stage-1: Complex FIR
- [ ] Stage-2: Complex IIR + SOS
- [ ] Stage-3: HackRF int8 FIR/IIR parity
- [ ] Stage-4: USRP B210 int16 LE FIR/IIR
- [ ] Stage-5: Коэффициенты и билдеры IQ-линейки
- [ ] Stage-6: Полная тестовая матрица и регрессия

## Acceptance (Epic)
- [ ] Доступны рабочие фильтры для complex/int8/int16
- [ ] Единый подход к расчёту коэффициентов сохранён
- [ ] CI содержит стабильную регрессию по IQ-линейке
- [ ] Есть пользовательская документация и примеры

## Технические риски
- Расхождение между complex/int ветками
- Устойчивость IIR высоких порядков
- Форматные ошибки int16 LE

## Митигации
- Golden vectors + cross-check pipeline
- SOS как обязательный путь для high-order
- Отдельные тесты endian/format
