# Низкочастотный эллиптический фильтр (порт C# EllipticLowPass)
import numpy as np
from math import ceil, sqrt, pi, tan, atan

from .elliptic_filter import EllipticFilter
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter


class EllipticLowPass(EllipticFilter):
    """
    Цифровой низкочастотный эллиптический фильтр (ФНЧ).

    Обеспечивает минимальный порядок при заданных ограничениях на пульсации
    в полосе пропускания и ослабление в полосе заграждения.
    """

    @staticmethod
    def get_order(dt: float, fp: float, fs: float,
                  Gp: float = 0.891250938, Gs: float = 0.01) -> int:
        """Минимальный порядок эллиптического ФНЧ по формуле через K и K'."""
        k_eps = sqrt((1.0 / (Gp * Gp) - 1.0) / (1.0 / (Gs * Gs) - 1.0))
        pi_dt = pi * dt
        kW = tan(fs * pi_dt) / tan(fp * pi_dt)
        K = EllipticFilter.K
        T = EllipticFilter.T
        return int(ceil(T(k_eps) * K(kW) / (K(k_eps) * T(kW))))

    @staticmethod
    def _initialize(spec: Specification):
        """
        Коэффициенты ФНЧ эллиптического фильтра.

        Алгоритм
        --------
        1. N = ceil(T(1/kEps)·K(1/kW) / (K(1/kEps)·T(1/kW))).
        2. Нули и полюса нормированного эллиптического прототипа.
        3. Масштабируем на Wp и преобразуем билинейным методом.
        4. Нечётный порядок: добавляем нуль z=-1.
        5. Нормировка по DC: k0/|prod(1-z_zero)/prod(1-z_pole)|.
        """
        K = EllipticFilter.K
        T = EllipticFilter.T

        # обратные модули (LP-прототип: kW и kEps перевёрнуты относительно BS)
        k_W = 1.0 / spec.kW
        k_eps = 1.0 / spec.kEps

        N = int(ceil(T(k_eps) * K(k_W) / (K(k_eps) * T(k_W))))

        zeros, poles = EllipticFilter.get_normed_zeros_poles(N, spec.EpsP, spec.EpsS)

        z_poles = DigitalFilter.to_z_array(poles, spec.dt, spec.Wp)

        # z-нули из аналоговых нулей; нечётный порядок → добавляем z=-1
        z_zeros_list = list(DigitalFilter.to_z_array(zeros, spec.dt, spec.Wp))
        if N % 2 != 0:
            z_zeros_list = [-1.0 + 0j] + z_zeros_list
        z_zeros = np.array(z_zeros_list)

        # нормировка усиления на нулевой частоте (z=1)
        k_z = np.prod([1.0 - z for z in z_zeros])
        k_p = np.prod([1.0 - z for z in z_poles])
        k0 = 1.0 if N % 2 != 0 else (1.0 / sqrt(1.0 + spec.EpsP ** 2))
        g_norm = k0 * abs(k_p / k_z)

        B = np.poly(z_zeros).real * g_norm
        A = np.poly(z_poles).real
        return A, B

    def __init__(self, dt: float, fp: float, fs: float,
                 Gp: float = 0.891250938, Gs: float = 0.01):
        """
        Создаёт ФНЧ эллиптический.

        Параметры
        ---------
        dt : период дискретизации (с)
        fp : граничная частота полосы пропускания (Гц)
        fs : граничная частота полосы заграждения (Гц)
        Gp : коэффициент передачи на границе пропускания
        Gs : коэффициент передачи на границе заграждения
        """
        spec = Specification(dt, fp, fs, Gp, Gs)
        A, B = EllipticLowPass._initialize(spec)
        super().__init__(B, A, spec)
