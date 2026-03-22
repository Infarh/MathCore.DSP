# Низкочастотный фильтр Баттерворта (порт C# ButterworthLowPass)
import numpy as np
from math import ceil, log, tan, atan, pi

from .butterworth_filter import ButterworthFilter
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter


class ButterworthLowPass(ButterworthFilter):
    """
    Цифровой низкочастотный фильтр Баттерворта (ФНЧ).

    Проектируется через аналоговый ФНЧ-прототип с последующим
    билинейным z-преобразованием.
    """

    # ── Вспомогательные статические методы ───────────────────────────────────

    @staticmethod
    def get_order(dt: float, fp: float, fs: float,
                  Gp: float = 0.891250938, Gs: float = 0.01) -> int:
        """Минимальный порядок фильтра по заданной спецификации."""
        k_eps2 = (1.0 / (Gs * Gs) - 1.0) / (1.0 / (Gp * Gp) - 1.0)
        pi_dt = pi * dt
        kW = tan(fs * pi_dt) / tan(fp * pi_dt)
        return int(ceil(0.5 * log(k_eps2) / log(kW)))

    @staticmethod
    def get_frequency_stop(dt: float, fp: float, Order: int,
                           Gp: float = 0.891250938, Gs: float = 0.01) -> float:
        """Частота заграждения для заданного порядка фильтра."""
        k_eps2 = (1.0 / (Gs * Gs) - 1.0) / (1.0 / (Gp * Gp) - 1.0)
        pi_dt = pi * dt
        return atan(k_eps2 ** (0.5 / Order) * tan(pi_dt * fp)) / pi_dt

    @staticmethod
    def get_frequency_pass(dt: float, Order: int, fs: float,
                           Gp: float = 0.891250938, Gs: float = 0.01) -> float:
        """Частота пропускания для заданного порядка фильтра."""
        k_eps2 = (1.0 / (Gs * Gs) - 1.0) / (1.0 / (Gp * Gp) - 1.0)
        pi_dt = pi * dt
        return atan(k_eps2 ** (-0.5 / Order) * tan(pi_dt * fs)) / pi_dt

    @staticmethod
    def get_poles(Order: int, Gp: float = 0.891250938) -> list:
        """Нормированные полюса ФНЧ Баттерворта заданного порядка."""
        return ButterworthFilter.get_norm_poles_gp(Order, Gp)

    # ── Расчёт коэффициентов ──────────────────────────────────────────────────

    @staticmethod
    def get_polynoms(spec: Specification):
        """
        Коэффициенты передаточной функции H(z) = B(z)/A(z).

        Алгоритм
        --------
        1. Вычисляем порядок N = ceil(log(kEps) / log(kW)).
        2. Получаем нормированные полюса Баттерворта.
        3. Масштабируем полюса на Wp (LP→LP трансформация).
        4. Билинейное z-преобразование: p → z.
        5. Нормирующий коэффициент g = prod((1-zi)/2).real.
        6. Числитель B — биномиальные коэффициенты * g.
        7. Знаменатель A = np.poly(z_poles).real.
        """
        N = int(ceil(log(spec.kEps) / log(spec.kW)))

        poles = ButterworthFilter.get_norm_poles(N, spec.EpsP)

        # LP→LP трансформация + билинейное z-преобразование
        z_poles = DigitalFilter.to_z_array(poles, spec.dt, spec.Wp)

        # нормирующий коэффициент: усиление на нулевой частоте (z=1)
        g_norm = np.real(np.prod([(1.0 - z) / 2.0 for z in z_poles]))

        # числитель: биномиальные коэффициенты C(N,k) * g_norm
        B = g_norm * AnalogBasedFilter.get_binomial_coefficients(N)

        # знаменатель: полином из корней z_poles
        A = np.poly(z_poles).real

        return A, B

    # ── Конструктор ───────────────────────────────────────────────────────────

    def __init__(self, dt: float, fp: float, fs: float,
                 Gp: float = 0.891250938, Gs: float = 0.01):
        """
        Создаёт ФНЧ Баттерворта.

        Параметры
        ---------
        dt : период дискретизации (с)
        fp : граничная частота полосы пропускания (Гц)
        fs : граничная частота полосы заграждения (Гц)
        Gp : коэффициент передачи на границе пропускания (-1 дБ по умолчанию)
        Gs : коэффициент передачи на границе заграждения (-40 дБ по умолчанию)
        """
        spec = Specification(dt, fp, fs, Gp, Gs)
        A, B = ButterworthLowPass.get_polynoms(spec)
        super().__init__(B, A, spec)
