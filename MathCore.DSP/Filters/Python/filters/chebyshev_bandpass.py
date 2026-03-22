# Полосовой фильтр Чебышева (порт C# ChebyshevBandPass)
import numpy as np
from math import ceil, pi, fabs

from .chebyshev_filter import ChebyshevFilter, ChebyshevType
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter

pi2 = 2.0 * pi


class ChebyshevBandPass(ChebyshevFilter):
    """
    Цифровой полосовой фильтр Чебышева (ППФ).

    Нули в z-области: N штук z=+1 (ФНЧ-прототип тип I)
    или из нулей Чебышева II рода (тип II).
    """

    @staticmethod
    def _get_specification(dt: float, fsl: float, fpl: float, fph: float, fsh: float,
                           Gp: float, Gs: float) -> Specification:
        """Спецификация ППФ через эквивалентный ФНЧ-прототип."""
        Fsl = DigitalFilter.to_digital_frequency(fsl, dt)
        Fpl = DigitalFilter.to_digital_frequency(fpl, dt)
        Fph = DigitalFilter.to_digital_frequency(fph, dt)
        Fsh = DigitalFilter.to_digital_frequency(fsh, dt)

        Wsl, Wpl = pi2 * Fsl, pi2 * Fpl
        Wph, Wsh = pi2 * Fph, pi2 * Fsh

        Wc = Wpl * Wph
        dW = Wph - Wpl

        Ws = Wsh if (Wc / Wsh > Wsl) else Wsl
        # W1 = |(Wc - Ws²) / (dW·Ws)| — прототип ФНЧ с fp=1/pi2, fs=W1/pi2
        W1 = fabs((Wc - Ws ** 2) / (dW * Ws))
        Fp = 1.0 / pi2
        Fs = W1 / pi2

        fp = DigitalFilter.to_analog_frequency(Fp, dt)
        fs = DigitalFilter.to_analog_frequency(Fs, dt)
        return Specification(dt, fp, fs, Gp, Gs)

    @staticmethod
    def _center_z(fpl: float, fph: float, dt: float) -> complex:
        """Значение z на центральной частоте полосы пропускания."""
        Fpl = DigitalFilter.to_digital_frequency(fpl, dt)
        Fph = DigitalFilter.to_digital_frequency(fph, dt)
        Fp0 = np.sqrt(Fpl * Fph)
        ffp0 = DigitalFilter.to_analog_frequency(Fp0, dt)
        return np.exp(1j * pi2 * ffp0 * dt)

    @staticmethod
    def _initialize_i(fpl: float, fph: float, spec: Specification):
        """Коэффициенты ППФ Чебышева I рода."""
        dt = spec.dt
        Wpl = pi2 * DigitalFilter.to_digital_frequency(fpl, dt)
        Wph = pi2 * DigitalFilter.to_digital_frequency(fph, dt)

        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        poles = ChebyshevFilter.get_normed_poles_i(N, spec.EpsP)

        # LP→BP трансформация полюсов
        ppf_poles = AnalogBasedFilter.transform_to_band_pass_w(poles, Wpl, Wph)

        # нули ФНЧ I рода в BP → z=+1 и z=-1 (N штук каждого)
        z_zeros = np.array([1.0 + 0j] * N + [-1.0 + 0j] * N)
        z_poles = DigitalFilter.to_z_array(ppf_poles, dt)

        # нормировка по центральной частоте
        z0 = ChebyshevBandPass._center_z(fpl, fph, dt)
        norm_0 = ((z0 - 1.0) * (z0 + 1.0)) ** N
        norm_p = np.prod([z0 - z for z in z_poles])

        k0 = spec.Gp if N % 2 == 0 else abs(z0)
        g_norm = k0 * abs(norm_p / norm_0)

        B = np.poly(z_zeros).real * g_norm
        A = np.poly(z_poles).real
        return A, B

    @staticmethod
    def _initialize_ii(fpl: float, fph: float, spec: Specification):
        """Коэффициенты ППФ Чебышева II рода."""
        dt = spec.dt
        Wpl = pi2 * DigitalFilter.to_digital_frequency(fpl, dt)
        Wph = pi2 * DigitalFilter.to_digital_frequency(fph, dt)

        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        zeros, poles = ChebyshevFilter.get_normed_poles_ii(N, spec.EpsS)

        ppf_zeros = AnalogBasedFilter.transform_to_band_pass_w(zeros, Wpl, Wph)
        if N % 2 != 0:
            ppf_zeros = list(ppf_zeros) + [0.0 + 0j]  # дополнительный нуль
        ppf_poles = AnalogBasedFilter.transform_to_band_pass_w(poles, Wpl, Wph)

        z_zeros_list = list(DigitalFilter.to_z_array(ppf_zeros, dt))
        if N % 2 != 0:
            z_zeros_list.append(-1.0 + 0j)
        z_zeros = np.array(z_zeros_list)
        z_poles = DigitalFilter.to_z_array(ppf_poles, dt)

        z0 = ChebyshevBandPass._center_z(fpl, fph, dt)
        norm_0 = np.prod([z0 - z for z in z_zeros])
        norm_p = np.prod([z0 - z for z in z_poles])
        g_norm = (spec.Gp if N % 2 == 0 else 1.0) * (norm_p / norm_0).real

        B = np.poly(z_zeros).real * g_norm
        A = np.poly(z_poles).real
        return A, B

    @staticmethod
    def _initialize_ii_corrected(fpl: float, fph: float, spec: Specification):
        """Коэффициенты ППФ Чебышева II рода с коррекцией частот."""
        dt = spec.dt
        Wpl = pi2 * DigitalFilter.to_digital_frequency(fpl, dt)
        Wph = pi2 * DigitalFilter.to_digital_frequency(fph, dt)

        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        # W0 = kW для коррекции прототипа
        zeros, poles = ChebyshevFilter.get_normed_poles_ii(N, spec.EpsS, spec.kW)

        ppf_zeros = AnalogBasedFilter.transform_to_band_pass_w(zeros, Wpl, Wph)
        if N % 2 != 0:
            ppf_zeros = list(ppf_zeros) + [0.0 + 0j]
        ppf_poles = list(AnalogBasedFilter.transform_to_band_pass_w(poles, Wpl, Wph))

        z_zeros_list = list(DigitalFilter.to_z_array(ppf_zeros, dt))
        if N % 2 != 0:
            z_zeros_list.append(-1.0 + 0j)
        z_zeros = np.array(z_zeros_list)
        z_poles = DigitalFilter.to_z_array(ppf_poles, dt)

        z0 = ChebyshevBandPass._center_z(fpl, fph, dt)
        norm_0 = np.prod([z0 - z for z in z_zeros])
        norm_p = np.prod([z0 - z for z in z_poles])
        k0 = spec.Gp if N % 2 == 0 else abs(z0)
        g_norm = k0 * abs(norm_p / norm_0)

        B = np.poly(z_zeros).real * g_norm
        A = np.poly(z_poles).real
        return A, B

    def __init__(self, dt: float, fsl: float, fpl: float, fph: float, fsh: float,
                 Gp: float = 0.891250938, Gs: float = 0.01,
                 filter_type: ChebyshevType = ChebyshevType.I):
        """
        Создаёт ППФ Чебышева.

        Параметры
        ---------
        dt          : период дискретизации (с)
        fsl         : нижняя граница полосы заграждения (Гц)
        fpl         : нижняя граница полосы пропускания (Гц)
        fph         : верхняя граница полосы пропускания (Гц)
        fsh         : верхняя граница полосы заграждения (Гц)
        Gp          : коэффициент передачи на границе пропускания
        Gs          : коэффициент передачи на границе заграждения
        filter_type : тип фильтра
        """
        spec = ChebyshevBandPass._get_specification(dt, fsl, fpl, fph, fsh, Gp, Gs)
        init = {
            ChebyshevType.I:           lambda s: ChebyshevBandPass._initialize_i(fpl, fph, s),
            ChebyshevType.II:          lambda s: ChebyshevBandPass._initialize_ii(fpl, fph, s),
            ChebyshevType.IICorrected: lambda s: ChebyshevBandPass._initialize_ii_corrected(fpl, fph, s),
        }
        A, B = init[filter_type](spec)
        super().__init__(B, A, spec, filter_type)
