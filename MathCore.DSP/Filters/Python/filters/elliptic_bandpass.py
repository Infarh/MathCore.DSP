# Полосовой эллиптический фильтр (порт C# EllipticBandPass)
import numpy as np
from math import ceil, pi

from .elliptic_filter import EllipticFilter
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter
from scipy.special import ellipk as scipy_ellipk

pi2 = 2.0 * pi


class EllipticBandPass(EllipticFilter):
    """
    Цифровой полосовой эллиптический фильтр (ППФ).

    Нормировка выполняется по центральной частоте полосы пропускания.
    """

    @staticmethod
    def _get_specification(dt: float, fsl: float, fpl: float, fph: float, fsh: float,
                           Gp: float, Gs: float) -> Specification:
        """Спецификация ППФ через эквивалентный ФНЧ-прототип."""
        from math import fabs
        Fsl = DigitalFilter.to_digital_frequency(fsl, dt)
        Fpl = DigitalFilter.to_digital_frequency(fpl, dt)
        Fph = DigitalFilter.to_digital_frequency(fph, dt)
        Fsh = DigitalFilter.to_digital_frequency(fsh, dt)

        Wsl, Wpl = pi2 * Fsl, pi2 * Fpl
        Wph, Wsh = pi2 * Fph, pi2 * Fsh

        Wc = Wpl * Wph
        dW = Wph - Wpl

        Ws = Wsh if (Wc / Wsh > Wsl) else Wsl
        Fp = fabs(dW * Ws / (Wc - Ws ** 2)) / pi2
        Fs = 1.0 / pi2

        fp = DigitalFilter.to_analog_frequency(Fp, dt)
        fs = DigitalFilter.to_analog_frequency(Fs, dt)
        return Specification(dt, fp, fs, Gp, Gs)

    @staticmethod
    def _initialize(fsl: float, fpl: float, fph: float, fsh: float, spec: Specification):
        """
        Коэффициенты ППФ эллиптического фильтра.

        Алгоритм
        --------
        1. N через эллиптические интегралы (с инвертированными модулями).
        2. Нули и полюса LP-прототипа.
        3. LP→BP трансформация.
        4. Нечётный порядок → нуль в z=-1.
        5. Нормировка по центральной частоте полосы пропускания.
        """
        dt = spec.dt
        Wpl = pi2 * DigitalFilter.to_digital_frequency(fpl, dt)
        Wph = pi2 * DigitalFilter.to_digital_frequency(fph, dt)

        k_eps = 1.0 / spec.kEps
        k_W = 1.0 / spec.kW

        K_w = float(scipy_ellipk(k_W ** 2))
        T_w = float(scipy_ellipk(1.0 - k_W ** 2))
        K_eps = float(scipy_ellipk(k_eps ** 2))
        T_eps = float(scipy_ellipk(1.0 - k_eps ** 2))

        N = int(ceil(T_eps * K_w / (K_eps * T_w)))

        zeros, poles = EllipticFilter.get_normed_zeros_poles(N, spec.EpsP, spec.EpsS)

        ppf_zeros = AnalogBasedFilter.transform_to_band_pass_w(zeros, Wpl, Wph)
        ppf_poles = AnalogBasedFilter.transform_to_band_pass_w(poles, Wpl, Wph)

        z_zeros_list = list(DigitalFilter.to_z_array(ppf_zeros, dt))
        if N % 2 != 0:
            z_zeros_list = [-1.0 + 0j] + z_zeros_list
        z_zeros = np.array(z_zeros_list)
        z_poles = DigitalFilter.to_z_array(ppf_poles, dt)

        # нормировка по центральной частоте полосы
        Fpl_d = DigitalFilter.to_digital_frequency(fpl, dt)
        Fph_d = DigitalFilter.to_digital_frequency(fph, dt)
        Fp0 = np.sqrt(Fpl_d * Fph_d)
        ffp0 = DigitalFilter.to_analog_frequency(Fp0, dt)
        z0 = np.exp(1j * pi2 * ffp0 * dt)

        norm_0 = np.prod([z0 - z for z in z_zeros])
        norm_p = np.prod([z0 - z for z in z_poles])
        g_norm = (spec.Gp if N % 2 == 0 else abs(z0)) * abs(norm_p / norm_0)

        B = np.poly(z_zeros).real * g_norm
        A = np.poly(z_poles).real
        return A, B

    def __init__(self, dt: float, fsl: float, fpl: float, fph: float, fsh: float,
                 Gp: float = 0.891250938, Gs: float = 0.01):
        """
        Создаёт ППФ эллиптический.

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
        spec = EllipticBandPass._get_specification(dt, fsl, fpl, fph, fsh, Gp, Gs)
        A, B = EllipticBandPass._initialize(fsl, fpl, fph, fsh, spec)
        super().__init__(B, A, spec)
