﻿<style>
.markdown-body h1{
    color:red;
}

.markdown-body h2{
    color:blue;
}
</style>

# Цифровые фильтры

Спад АЧХ расчитывается в дБ/октава (на удвоенное знечение частоты).

Фильтр <span style="color:red">**Проверить!!!**</span>.
- 1-порядка - спад 3дБ/окатва (на удвоенной частоте среда затухание -6 дБ)
- 2-порядка - сапд 6дБ/октава (на удвоенной частоте среза затухание -12дБ)
- 3 порядка - спад 12дБ/октава (-18 дБ на удвоенной частоте среза)

## КИХ-фильтры (FIR)

- Отдельно управляются АЧХ и ФЧХ

## БИХ-фильтры (IIR)

### Фильтр Баттерворта

- Самый простой и распространённый
- Гладкая АЧХ
- бОльший порядок для заданной полосы
- Сдвигает фазу на 45°/1-порядок (?)

### Эллиптический фильтр

### Фильтр Чебышева

#### I-рода

#### II-рода

- Нужен для гарантированного подавления полосы частот

#### II-рода коректированный

- ОБеспечивает точное попадание в требуемые частоту среза и частоту заграждения

### Фильтр Бесселя


