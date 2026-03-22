# Высокочастотный эллиптический фильтр (порт C# EllipticHighPass)
import numpy as np
from math import ceil

from .elliptic_filter import EllipticFilter
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter


class EllipticHighPass(EllipticFilter):
    """
    Цифровой высокочастотный эллиптический фильтр (ФВЧ).

    LP→HP трансформация нормированных нулей и полюсов: p → Ws/p.
    """

    @staticmethod
    def get_polynoms(spec: Specification):
        """
        Коэффициенты ФВЧ эллиптического фильтра.

        Алгоритм
        --------
        1. N = ceil(T(1/kEps)·K(1/kw) / (K(1/kEps)·T(1/kw))).
        2. Нули и полюса нормированного прототипа.
        3. LP→HP трансформация: p → Ws/p.
        4. Нечётный порядок → дополнительный нуль в 0.
        5. Нормировка по Найквисту: prod((1+zi)/2).
        """
        K = EllipticFilter.K
        T = EllipticFilter.T

        k_w = 1.0 / spec.kw
        k_eps = 1.0 / spec.kEps

        N = int(ceil(T(k_eps) * K(k_w) / (K(k_eps) * T(k_w))))

        zeros, poles = EllipticFilter.get_normed_zeros_poles(N, spec.EpsP, spec.EpsS)

        # LP→HP трансформация
        hp_zeros = AnalogBasedFilter.transform_to_high_pass_w(zeros, spec.Ws)
        if N % 2 != 0:
            hp_zeros = [0.0 + 0j] + list(hp_zeros)  # нуль на нулевой частоте
        hp_poles = AnalogBasedFilter.transform_to_high_pass_w(poles, spec.Ws)

        z_zeros = DigitalFilter.to_z_array(hp_zeros, spec.dt)
        z_poles = DigitalFilter.to_z_array(hp_poles, spec.dt)

        # нормировка по Найквисту (z = -1)
        k_zeros = np.prod([1.0 + z for z in z_zeros]).real
        k_poles = np.prod([1.0 + z for z in z_poles]).real
        g_norm = (spec.Gp if N % 2 == 0 else 1.0) * k_poles / k_zeros

        B = np.poly(z_zeros).real * g_norm
        A = np.poly(z_poles).real
        return A, B

    def __init__(self, dt: float, fs: float, fp: float,
                 Gp: float = 0.891250938, Gs: float = 0.01):
        """
        Создаёт ФВЧ эллиптический.

        Параметры
        ---------
        dt : период дискретизации (с)
        fs : частота заграждения (Гц), fs < fp
        fp : граничная частота полосы пропускания (Гц)
        Gp : коэффициент передачи на границе пропускания
        Gs : коэффициент передачи на границе заграждения
        """
        spec = AnalogBasedFilter.get_specification(dt, fs, fp, Gp, Gs)
        A, B = EllipticHighPass.get_polynoms(spec)
        super().__init__(B, A, spec)
