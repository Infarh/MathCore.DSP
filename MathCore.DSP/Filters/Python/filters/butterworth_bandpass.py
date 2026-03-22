# Полосовой фильтр Баттерворта (порт C# ButterworthBandPass)
import numpy as np
from math import ceil, log, pi, tan, atan, sqrt, fabs

from .butterworth_filter import ButterworthFilter
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter

pi2 = 2.0 * pi


class ButterworthBandPass(ButterworthFilter):
    """
    Цифровой полосовой фильтр Баттерворта (ППФ).

    Нули в z-области: N нулей в точке z=+1 (постоянная составляющая)
    и N нулей в точке z=-1 (частота Найквиста).
    """

    @staticmethod
    def _get_specification(dt: float, fsl: float, fpl: float, fph: float, fsh: float,
                           Gp: float, Gs: float) -> Specification:
        """Формирует спецификацию ППФ через эквивалентный ФНЧ-прототип."""
        Fsl = DigitalFilter.to_digital_frequency(fsl, dt)
        Fpl = DigitalFilter.to_digital_frequency(fpl, dt)
        Fph = DigitalFilter.to_digital_frequency(fph, dt)
        Fsh = DigitalFilter.to_digital_frequency(fsh, dt)

        Wsl, Wpl = pi2 * Fsl, pi2 * Fpl
        Wph, Wsh = pi2 * Fph, pi2 * Fsh

        Wc = Wpl * Wph   # квадрат центральной частоты
        dW = Wph - Wpl   # ширина полосы пропускания

        # выбираем опорную частоту: верхнюю или нижнюю границу заграждения
        Ws = Wsh if (Wc / Wsh > Wsl) else Wsl
        Fp = fabs(dW * Ws / (Wc - Ws ** 2)) / pi2
        Fs = 1.0 / pi2  # нормировка прототипа

        fp = DigitalFilter.to_analog_frequency(Fp, dt)
        fs = DigitalFilter.to_analog_frequency(Fs, dt)
        return Specification(dt, fp, fs, Gp, Gs)

    @staticmethod
    def _initialize(fpl: float, fph: float, spec: Specification):
        """
        Расчёт коэффициентов ППФ Баттерворта.

        Алгоритм
        --------
        1. N = ceil(log(kEps)/log(kW)).
        2. Нормированные полюса ФНЧ-прототипа.
        3. LP→BP трансформация: каждый полюс даёт 2 полюса ППФ.
        4. Нули в z-области: N×(z=+1) + N×(z=-1).
        5. Нормировка усиления по аналоговым ВЧ-полюсам:
           g = ∏ fd2·dW / ((fd2-p1)(fd2-p2)) / EpsP.
        """
        dt = spec.dt
        Wpl = pi2 * DigitalFilter.to_digital_frequency(fpl, dt)
        Wph = pi2 * DigitalFilter.to_digital_frequency(fph, dt)
        dW = Wph - Wpl

        N = int(ceil(log(spec.kEps) / log(spec.kW)))

        poles = ButterworthFilter.get_norm_poles(N, spec.EpsP)

        # LP→BP трансформация (каждый LP-полюс → 2 BP-полюса)
        ppf_poles = AnalogBasedFilter.transform_to_band_pass_w(poles, Wpl, Wph)

        # нули в z-области: N штук z=+1 и N штук z=-1
        z_zeros = np.array([1.0 + 0j] * N + [-1.0 + 0j] * N)

        # билинейное z-преобразование аналоговых BP-полюсов
        z_poles = DigitalFilter.to_z_array(ppf_poles, dt)

        # нормировка усиления по парам аналоговых BP-полюсов
        fd2 = 2.0 * spec.fd
        g_norm = complex(1.0, 0.0)
        for i in range(N):
            p1, p2 = ppf_poles[2 * i], ppf_poles[2 * i + 1]
            g_norm *= fd2 * dW / ((fd2 - p1) * (fd2 - p2))
        G_norm = g_norm / spec.EpsP

        B = np.poly(z_zeros).real * G_norm.real
        A = np.poly(z_poles).real

        return A, B

    def __init__(self, dt: float, fsl: float, fpl: float, fph: float, fsh: float,
                 Gp: float = 0.891250938, Gs: float = 0.01):
        """
        Создаёт ППФ Баттерворта.

        Параметры
        ---------
        dt  : период дискретизации (с)
        fsl : нижняя граница полосы заграждения (Гц)
        fpl : нижняя граница полосы пропускания (Гц)
        fph : верхняя граница полосы пропускания (Гц)
        fsh : верхняя граница полосы заграждения (Гц)
        Gp  : коэффициент передачи на границе пропускания
        Gs  : коэффициент передачи на границе заграждения
        """
        spec = ButterworthBandPass._get_specification(dt, fsl, fpl, fph, fsh, Gp, Gs)
        A, B = ButterworthBandPass._initialize(fpl, fph, spec)
        super().__init__(B, A, spec)
