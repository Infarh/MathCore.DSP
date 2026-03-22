# Полосозаграждающий фильтр Баттерворта (порт C# ButterworthBandStop)
import numpy as np
from math import ceil, log, pi, fabs

from .butterworth_filter import ButterworthFilter
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter

pi2 = 2.0 * pi


class ButterworthBandStop(ButterworthFilter):
    """
    Цифровой полосозаграждающий фильтр Баттерворта (ПЗФ).

    Нули в аналоговой области: 2N штук на мнимой оси в точке ±j·√Wc,
    что соответствует центру полосы заграждения.
    """

    @staticmethod
    def _get_specification(dt: float, fpl: float, fsl: float, fsh: float, fph: float,
                           Gp: float, Gs: float) -> Specification:
        """Формирует спецификацию ПЗФ через эквивалентный ФНЧ-прототип."""
        Fpl = DigitalFilter.to_digital_frequency(fpl, dt)
        Fsl = DigitalFilter.to_digital_frequency(fsl, dt)
        Fsh = DigitalFilter.to_digital_frequency(fsh, dt)
        Fph = DigitalFilter.to_digital_frequency(fph, dt)

        Wpl, Wsl = pi2 * Fpl, pi2 * Fsl
        Wsh, Wph = pi2 * Fsh, pi2 * Fph

        Wc = Wsl * Wsh   # квадрат центральной частоты полосы заграждения
        dW = Wsh - Wsl   # ширина полосы заграждения

        # опорная частота: верхняя или нижняя граница пропускания
        Wp = Wph if (Wc / Wph > Wpl) else Wpl
        Fp = fabs(dW * Wp / (Wc - Wp ** 2)) / pi2
        Fs = 1.0 / pi2

        fp = DigitalFilter.to_analog_frequency(Fp, dt)
        fs = DigitalFilter.to_analog_frequency(Fs, dt)
        return Specification(dt, fp, fs, Gp, Gs)

    @staticmethod
    def _initialize(fpl: float, fsl: float, fsh: float, fph: float, spec: Specification):
        """
        Расчёт коэффициентов ПЗФ Баттерворта.

        Алгоритм
        --------
        1. Вычисляем опорный коэффициент W0 и порядок N.
        2. Нормированные полюса ФНЧ-прототипа с масштабом W0.
        3. Нули в аналоговой области: 2N штук ±j·√Wc.
        4. LP→BS трансформация полюсов.
        5. Билинейное z-преобразование нулей и полюсов.
        6. Нормировка на нулевой частоте (z=1).
        """
        dt = spec.dt
        Wpl = pi2 * DigitalFilter.to_digital_frequency(fpl, dt)
        Wsl = pi2 * DigitalFilter.to_digital_frequency(fsl, dt)
        Wsh = pi2 * DigitalFilter.to_digital_frequency(fsh, dt)
        Wph = pi2 * DigitalFilter.to_digital_frequency(fph, dt)

        Wc = Wsl * Wsh
        dW = Wsh - Wsl
        sqrt_wc = np.sqrt(Wc)

        # опорная частота и масштаб для ФНЧ-прототипа
        Wp = Wph if (Wc / Wph > Wpl) else Wpl
        W0 = fabs(dW * Wp / (Wc - Wp ** 2))

        N = int(ceil(log(spec.kEps) / log(spec.kW)))

        poles = ButterworthFilter.get_norm_poles(N, spec.EpsP, W0)

        # нули в аналоговой области: N пар ±j·√Wc
        pzf_zeros = [complex(0.0, sqrt_wc if i % 2 == 0 else -sqrt_wc)
                     for i in range(2 * N)]

        # LP→BS трансформация полюсов
        pzf_poles = AnalogBasedFilter.transform_to_band_stop_w(poles, Wsl, Wsh)

        # билинейное z-преобразование
        z_zeros = DigitalFilter.to_z_array(pzf_zeros, dt)
        z_poles = DigitalFilter.to_z_array(pzf_poles, dt)

        # нормировка: усиление на нулевой частоте (z=1)
        gain_zeros = np.prod([1.0 - z for z in z_zeros])
        gain_poles = np.prod([1.0 - z for z in z_poles])
        k0 = 1.0 if N % 2 != 0 else spec.Gp
        g_norm = k0 / abs(gain_zeros / gain_poles)

        B = np.poly(z_zeros).real * g_norm
        A = np.poly(z_poles).real

        return A, B

    def __init__(self, dt: float, fpl: float, fsl: float, fsh: float, fph: float,
                 Gp: float = 0.891250938, Gs: float = 0.01):
        """
        Создаёт ПЗФ Баттерворта.

        Параметры
        ---------
        dt  : период дискретизации (с)
        fpl : нижняя граница полосы пропускания (Гц)
        fsl : нижняя граница полосы заграждения (Гц)
        fsh : верхняя граница полосы заграждения (Гц)
        fph : верхняя граница полосы пропускания (Гц)
        Gp  : коэффициент передачи на границе пропускания
        Gs  : коэффициент передачи на границе заграждения
        """
        spec = ButterworthBandStop._get_specification(dt, fpl, fsl, fsh, fph, Gp, Gs)
        A, B = ButterworthBandStop._initialize(fpl, fsl, fsh, fph, spec)
        super().__init__(B, A, spec)
