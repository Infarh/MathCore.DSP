# Высокочастотный фильтр Чебышева (порт C# ChebyshevHighPass)
import numpy as np
from math import ceil

from .chebyshev_filter import ChebyshevFilter, ChebyshevType
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter


class ChebyshevHighPass(ChebyshevFilter):
    """
    Цифровой высокочастотный фильтр Чебышева (ФВЧ).

    Для LP→HP трансформации используются нормированные полюса ФНЧ Чебышева
    и трансформация p → Ws/p, где Ws — циклическая частота среза.
    """

    @staticmethod
    def _initialize_i(spec: Specification):
        """Коэффициенты ФВЧ Чебышева I рода."""
        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        poles = ChebyshevFilter.get_normed_poles_i(N, spec.EpsP)

        # LP→HP трансформация: p → Ws/p
        hp_poles = AnalogBasedFilter.transform_to_high_pass_w(poles, spec.Ws)

        z_poles = DigitalFilter.to_z_array(hp_poles, spec.dt)

        # нормировка по частоте Найквиста (z = -1)
        k_poles = np.prod([(1.0 + z) / 2.0 for z in z_poles]).real
        g_norm = (spec.Gp if N % 2 == 0 else 1.0) * k_poles

        binomial = AnalogBasedFilter.get_binomial_coefficients(N)
        B = np.array([b * (g_norm if i % 2 == 0 else -g_norm)
                      for i, b in enumerate(binomial)])
        A = np.poly(z_poles).real

        return A, B

    @staticmethod
    def _initialize_ii(spec: Specification):
        """Коэффициенты ФВЧ Чебышева II рода."""
        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        zeros, poles = ChebyshevFilter.get_normed_poles_ii(N, spec.EpsS)

        # LP→HP трансформация нулей и полюсов
        hp_zeros = AnalogBasedFilter.transform_to_high_pass_w(zeros, spec.Ws)
        if N % 2 != 0:
            hp_zeros = [0.0 + 0j] + list(hp_zeros)  # дополнительный нуль в 0
        hp_poles = AnalogBasedFilter.transform_to_high_pass_w(poles, spec.Ws)

        z_zeros = DigitalFilter.to_z_array(hp_zeros, spec.dt)
        z_poles = DigitalFilter.to_z_array(hp_poles, spec.dt)

        B = np.poly(z_zeros).real
        A = np.poly(z_poles).real

        # нормировка по частоте Найквиста
        k_zeros = np.prod([1.0 + z for z in z_zeros]).real
        k_poles = np.prod([1.0 + z for z in z_poles]).real
        g_norm = k_poles / k_zeros

        B = B * g_norm
        return A, B

    @staticmethod
    def _initialize_ii_corrected(spec: Specification):
        """Коэффициенты ФВЧ Чебышева II рода с коррекцией частот."""
        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        # W0 = kw для коррекции прототипа
        zeros, poles = ChebyshevFilter.get_normed_poles_ii(N, spec.EpsS, spec.kw)

        hp_zeros = AnalogBasedFilter.transform_to_high_pass_w(zeros, spec.Ws)
        if N % 2 != 0:
            hp_zeros = [0.0 + 0j] + list(hp_zeros)
        hp_poles = AnalogBasedFilter.transform_to_high_pass_w(poles, spec.Ws)

        z_zeros = DigitalFilter.to_z_array(hp_zeros, spec.dt)
        z_poles = DigitalFilter.to_z_array(hp_poles, spec.dt)

        B = np.poly(z_zeros).real
        A = np.poly(z_poles).real

        k_zeros = np.prod([1.0 + z for z in z_zeros]).real
        k_poles = np.prod([1.0 + z for z in z_poles]).real
        g_norm = k_poles / k_zeros

        B = B * g_norm
        return A, B

    def __init__(self, dt: float, fs: float, fp: float,
                 Gp: float = 0.891250938, Gs: float = 0.01,
                 filter_type: ChebyshevType = ChebyshevType.I):
        """
        Создаёт ФВЧ Чебышева.

        Параметры
        ---------
        dt          : период дискретизации (с)
        fs          : частота заграждения (Гц), fs < fp
        fp          : граничная частота полосы пропускания (Гц)
        Gp          : коэффициент передачи на границе пропускания
        Gs          : коэффициент передачи на границе заграждения
        filter_type : тип фильтра
        """
        spec = AnalogBasedFilter.get_specification(dt, fs, fp, Gp, Gs)
        init = {
            ChebyshevType.I:           ChebyshevHighPass._initialize_i,
            ChebyshevType.II:          ChebyshevHighPass._initialize_ii,
            ChebyshevType.IICorrected: ChebyshevHighPass._initialize_ii_corrected,
        }
        A, B = init[filter_type](spec)
        super().__init__(B, A, spec, filter_type)
