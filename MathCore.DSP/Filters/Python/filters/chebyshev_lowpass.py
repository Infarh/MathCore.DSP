# Низкочастотный фильтр Чебышева (порт C# ChebyshevLowPass)
import numpy as np
from math import ceil

from .chebyshev_filter import ChebyshevFilter, ChebyshevType
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter


class ChebyshevLowPass(ChebyshevFilter):
    """
    Цифровой низкочастотный фильтр Чебышева (ФНЧ).

    Поддерживает три варианта:
    - Type I   : равноволновая пульсация в полосе пропускания
    - Type II  : равноволновая пульсация в полосе заграждения
    - IICorrected : II рода с точным попаданием в граничные частоты
    """

    # ── Инициализация для каждого типа ───────────────────────────────────────

    @staticmethod
    def _initialize_i(spec: Specification):
        """Коэффициенты ФНЧ Чебышева I рода."""
        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        poles = ChebyshevFilter.get_normed_poles_i(N, spec.EpsP)

        # LP-прототип → z-плоскость (с масштабом Wp)
        z_poles = DigitalFilter.to_z_array(poles, spec.dt, spec.Wp)

        A = np.poly(z_poles).real

        # нормировка: единичное усиление на нулевой частоте (чётный порядок → коррекция на Gp)
        k0 = 1.0 if N % 2 != 0 else spec.Gp
        prod_dc = np.prod([1.0 - z for z in z_poles]).real
        g_norm = k0 / (2.0 ** N / prod_dc)

        B = g_norm * AnalogBasedFilter.get_binomial_coefficients(N)
        return A, B

    @staticmethod
    def _initialize_ii(spec: Specification):
        """Коэффициенты ФНЧ Чебышева II рода."""
        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        zeros, poles = ChebyshevFilter.get_normed_poles_ii(N, spec.EpsS)

        z_zeros = DigitalFilter.to_z_array(zeros, spec.dt, spec.Wp)
        if N % 2 != 0:
            z_zeros = np.append([-1.0 + 0j], z_zeros)  # дополнительный нуль в z=-1
        z_poles = DigitalFilter.to_z_array(poles, spec.dt, spec.Wp)

        B = np.poly(z_zeros).real
        A = np.poly(z_poles).real

        # нормировка по постоянной составляющей
        k_zeros = np.prod([1.0 - z for z in z_zeros]).real
        k_poles = np.prod([1.0 - z for z in z_poles]).real
        g_norm = 1.0 / (k_zeros / k_poles)

        B = B * g_norm
        return A, B

    @staticmethod
    def _initialize_ii_corrected(spec: Specification):
        """Коэффициенты ФНЧ Чебышева II рода с коррекцией частот среза."""
        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        zeros, poles = ChebyshevFilter.get_normed_poles_ii(N, spec.EpsS)

        # масштабируем на Wp * kw для точного попадания в граничные частоты
        w_scale = spec.Wp * spec.kw
        z_zeros = DigitalFilter.to_z_array(zeros, spec.dt, w_scale)
        if N % 2 != 0:
            z_zeros = np.append([-1.0 + 0j], z_zeros)
        z_poles = DigitalFilter.to_z_array(poles, spec.dt, w_scale)

        B = np.poly(z_zeros).real
        A = np.poly(z_poles).real

        k_zeros = np.prod([1.0 - z for z in z_zeros]).real
        k_poles = np.prod([1.0 - z for z in z_poles]).real
        g_norm = 1.0 / (k_zeros / k_poles)

        B = B * g_norm
        return A, B

    # ── Конструктор ───────────────────────────────────────────────────────────

    def __init__(self, dt: float, fp: float, fs: float,
                 Gp: float = 0.891250938, Gs: float = 0.01,
                 filter_type: ChebyshevType = ChebyshevType.I):
        """
        Создаёт ФНЧ Чебышева.

        Параметры
        ---------
        dt          : период дискретизации (с)
        fp          : граничная частота полосы пропускания (Гц)
        fs          : граничная частота полосы заграждения (Гц)
        Gp          : коэффициент передачи на границе пропускания
        Gs          : коэффициент передачи на границе заграждения
        filter_type : тип фильтра (I, II или IICorrected)
        """
        spec = Specification(dt, fp, fs, Gp, Gs)
        init = {
            ChebyshevType.I:           ChebyshevLowPass._initialize_i,
            ChebyshevType.II:          ChebyshevLowPass._initialize_ii,
            ChebyshevType.IICorrected: ChebyshevLowPass._initialize_ii_corrected,
        }
        A, B = init[filter_type](spec)
        super().__init__(B, A, spec, filter_type)
