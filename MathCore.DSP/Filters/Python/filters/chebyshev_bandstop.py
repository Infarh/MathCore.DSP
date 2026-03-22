# Полосозаграждающий фильтр Чебышева (порт C# ChebyshevBandStop)
import numpy as np
from math import ceil, pi, fabs

from .chebyshev_filter import ChebyshevFilter, ChebyshevType
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter

pi2 = 2.0 * pi


class ChebyshevBandStop(ChebyshevFilter):
    """Цифровой полосозаграждающий фильтр Чебышева (ПЗФ)."""

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
    def _common_params(fpl: float, fsl: float, fsh: float, fph: float, dt: float):
        """Вычисляет общие параметры ПЗФ (Wc, dW, sqrt_wc, W0, Wp)."""
        Wpl = pi2 * DigitalFilter.to_digital_frequency(fpl, dt)
        Wsl = pi2 * DigitalFilter.to_digital_frequency(fsl, dt)
        Wsh = pi2 * DigitalFilter.to_digital_frequency(fsh, dt)
        Wph = pi2 * DigitalFilter.to_digital_frequency(fph, dt)

        Wc = Wsl * Wsh
        dW = Wsh - Wsl
        sqrt_wc = np.sqrt(Wc)

        Wp_ref = Wph if (Wc / Wph > Wpl) else Wpl
        W0 = fabs(dW * Wp_ref / (Wc - Wp_ref ** 2))
        return Wsl, Wsh, Wc, dW, sqrt_wc, W0

    @staticmethod
    def _normalize_bs(z_zeros, z_poles, N, spec):
        """Нормировка ПЗФ по постоянной составляющей."""
        gain_zeros = np.prod([1.0 - z for z in z_zeros])
        gain_poles = np.prod([1.0 - z for z in z_poles])
        k0 = 1.0 if N % 2 != 0 else spec.Gp
        return k0 / abs(gain_zeros / gain_poles)

    @staticmethod
    def _initialize_i(fpl: float, fsl: float, fsh: float, fph: float, spec: Specification):
        """Коэффициенты ПЗФ Чебышева I рода."""
        dt = spec.dt
        Wsl, Wsh, Wc, dW, sqrt_wc, W0 = ChebyshevBandStop._common_params(
            fpl, fsl, fsh, fph, dt)

        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        poles = ChebyshevFilter.get_normed_poles_i(N, spec.EpsP, W0)

        pzf_zeros = [complex(0.0, sqrt_wc if i % 2 == 0 else -sqrt_wc)
                     for i in range(2 * N)]
        pzf_poles = AnalogBasedFilter.transform_to_band_stop_w(poles, Wsl, Wsh)

        z_zeros = DigitalFilter.to_z_array(pzf_zeros, dt)
        z_poles = DigitalFilter.to_z_array(pzf_poles, dt)

        g_norm = ChebyshevBandStop._normalize_bs(z_zeros, z_poles, N, spec)
        return np.poly(z_poles).real, np.poly(z_zeros).real * g_norm

    @staticmethod
    def _initialize_ii(fpl: float, fsl: float, fsh: float, fph: float, spec: Specification):
        """Коэффициенты ПЗФ Чебышева II рода."""
        dt = spec.dt
        Wsl, Wsh, Wc, dW, sqrt_wc, W0 = ChebyshevBandStop._common_params(
            fpl, fsl, fsh, fph, dt)

        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(spec.kW)))

        zeros, poles = ChebyshevFilter.get_normed_poles_ii(N, spec.EpsS, W0)

        pzf_poles = AnalogBasedFilter.transform_to_band_stop_w(poles, Wsl, Wsh)
        pzf_zeros = list(AnalogBasedFilter.transform_to_band_stop_w(zeros, Wsl, Wsh))
        if N % 2 != 0:
            pzf_zeros += [complex(0.0, sqrt_wc), complex(0.0, -sqrt_wc)]

        z_zeros = DigitalFilter.to_z_array(pzf_zeros, dt)
        z_poles = DigitalFilter.to_z_array(pzf_poles, dt)

        gain_zeros = np.prod([1.0 - z for z in z_zeros])
        gain_poles = np.prod([1.0 - z for z in z_poles])
        g_norm = 1.0 / abs(gain_zeros / gain_poles)

        return np.poly(z_poles).real, np.poly(z_zeros).real * g_norm

    @staticmethod
    def _initialize_ii_corrected(fpl: float, fsl: float, fsh: float, fph: float,
                                 spec: Specification):
        """Коэффициенты ПЗФ Чебышева II рода с коррекцией частот."""
        dt = spec.dt
        Wsl, Wsh, Wc, dW, sqrt_wc, W0 = ChebyshevBandStop._common_params(
            fpl, fsl, fsh, fph, dt)

        N = int(ceil(ChebyshevFilter.arcch(spec.kEps) / ChebyshevFilter.arcch(1.0 / W0)))

        zeros, poles = ChebyshevFilter.get_normed_poles_ii(N, spec.EpsS)

        pzf_poles = AnalogBasedFilter.transform_to_band_stop_w(poles, Wsl, Wsh)
        pzf_zeros = list(AnalogBasedFilter.transform_to_band_stop_w(zeros, Wsl, Wsh))
        if N % 2 != 0:
            pzf_zeros += [complex(0.0, sqrt_wc), complex(0.0, -sqrt_wc)]

        z_zeros = DigitalFilter.to_z_array(pzf_zeros, dt)
        z_poles = DigitalFilter.to_z_array(pzf_poles, dt)

        gain_zeros = np.prod([1.0 - z for z in z_zeros])
        gain_poles = np.prod([1.0 - z for z in z_poles])
        g_norm = 1.0 / (gain_zeros / gain_poles).real

        return np.poly(z_poles).real, np.poly(z_zeros).real * g_norm

    def __init__(self, dt: float, fpl: float, fsl: float, fsh: float, fph: float,
                 Gp: float = 0.891250938, Gs: float = 0.01,
                 filter_type: ChebyshevType = ChebyshevType.I):
        """
        Создаёт ПЗФ Чебышева.

        Параметры
        ---------
        dt          : период дискретизации (с)
        fpl         : нижняя граница полосы пропускания (Гц)
        fsl         : нижняя граница полосы заграждения (Гц)
        fsh         : верхняя граница полосы заграждения (Гц)
        fph         : верхняя граница полосы пропускания (Гц)
        Gp          : коэффициент передачи на границе пропускания
        Gs          : коэффициент передачи на границе заграждения
        filter_type : тип фильтра
        """
        spec = ChebyshevBandStop._get_specification(dt, fpl, fsl, fsh, fph, Gp, Gs)
        init = {
            ChebyshevType.I:
                lambda s: ChebyshevBandStop._initialize_i(fpl, fsl, fsh, fph, s),
            ChebyshevType.II:
                lambda s: ChebyshevBandStop._initialize_ii(fpl, fsl, fsh, fph, s),
            ChebyshevType.IICorrected:
                lambda s: ChebyshevBandStop._initialize_ii_corrected(fpl, fsl, fsh, fph, s),
        }
        A, B = init[filter_type](spec)
        super().__init__(B, A, spec, filter_type)
