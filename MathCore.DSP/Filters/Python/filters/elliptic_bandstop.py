# Полосозаграждающий эллиптический фильтр (порт C# EllipticBandStop)
import numpy as np
from math import ceil, pi, fabs, sqrt

from .elliptic_filter import EllipticFilter
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter

pi2 = 2.0 * pi


class EllipticBandStop(EllipticFilter):
    """
    Цифровой полосозаграждающий эллиптический фильтр (ПЗФ).

    Нули в аналоговой области расположены на мнимой оси вблизи центра полосы заграждения.
    Нечётный порядок добавляет дополнительную пару нулей на центральной частоте.
    """

    @staticmethod
    def _get_specification(dt: float, fpl: float, fsl: float, fsh: float, fph: float,
                           Gp: float, Gs: float) -> Specification:
        """Спецификация ПЗФ через эквивалентный ФНЧ-прототип."""
        Fpl = DigitalFilter.to_digital_frequency(fpl, dt)
        Fsl = DigitalFilter.to_digital_frequency(fsl, dt)
        Fsh = DigitalFilter.to_digital_frequency(fsh, dt)
        Fph = DigitalFilter.to_digital_frequency(fph, dt)

        Wpl, Wsl = pi2 * Fpl, pi2 * Fsl
        Wsh, Wph = pi2 * Fsh, pi2 * Fph

        Wc = Wsl * Wsh
        dW = Wsh - Wsl

        Wp = Wph if (Wc / Wph > Wpl) else Wpl
        Fp = fabs(dW * Wp / (Wc - Wp ** 2)) / pi2
        Fs = 1.0 / pi2

        fp = DigitalFilter.to_analog_frequency(Fp, dt)
        fs = DigitalFilter.to_analog_frequency(Fs, dt)
        return Specification(dt, fp, fs, Gp, Gs)

    @staticmethod
    def _initialize(fsl: float, fsh: float, spec: Specification):
        """
        Коэффициенты ПЗФ эллиптического фильтра.

        Алгоритм
        --------
        1. N = ceil(T(1/kEps)·K(1/kW) / (K(1/kEps)·T(1/kW))).
        2. Нули и полюса прototипа масштабируются на W0=Wp.
        3. LP→BS трансформация через цифровые частоты Fsl, Fsh.
        4. Нечётный порядок → дополнительная пара нулей ±j·2π·√(fsl·fsh).
        5. Нормировка по нулевой частоте (DC).

        Примечание: нечётный дополнительный нуль использует исходные fsl, fsh,
        а не предыскажённые Fsl, Fsh — соответствует оригинальному C#-коду.
        """
        K = EllipticFilter.K
        T = EllipticFilter.T

        k_W = 1.0 / spec.kW
        k_eps = 1.0 / spec.kEps

        N = int(ceil(T(k_eps) * K(k_W) / (K(k_eps) * T(k_W))))

        # нули и полюса LP-прototипа масштабированы на Wp
        zeros, poles = EllipticFilter.get_normed_zeros_poles(N, spec.EpsP, spec.EpsS, spec.Wp)

        dt = spec.dt
        Fsl = DigitalFilter.to_digital_frequency(fsl, dt)
        Fsh = DigitalFilter.to_digital_frequency(fsh, dt)

        # LP→BS трансформация через предыскажённые частоты Fsl, Fsh
        pzf_zeros = AnalogBasedFilter.transform_to_band_stop(zeros, Fsl, Fsh)
        pzf_poles = AnalogBasedFilter.transform_to_band_stop(poles, Fsl, Fsh)

        # нечётный порядок: дополнительные нули на центральной частоте полосы заграждения
        # используются исходные fsl, fsh (не Fsl, Fsh) — соответствует C#-коду
        if N % 2 != 0:
            wc_sqrt = pi2 * sqrt(fsl * fsh)
            pzf_zeros = list(pzf_zeros) + [complex(0.0, +wc_sqrt), complex(0.0, -wc_sqrt)]

        z_zeros = DigitalFilter.to_z_array(pzf_zeros, dt)
        z_poles = DigitalFilter.to_z_array(pzf_poles, dt)

        # нормировка по DC (z=1)
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
        Создаёт ПЗФ эллиптический.

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
        spec = EllipticBandStop._get_specification(dt, fpl, fsl, fsh, fph, Gp, Gs)
        A, B = EllipticBandStop._initialize(fsl, fsh, spec)
        super().__init__(B, A, spec)
