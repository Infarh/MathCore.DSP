# Базовый класс фильтра Баттерворта (порт C# ButterworthFilter)
import numpy as np
from .analog_based_filter import AnalogBasedFilter, Specification

_pi05 = np.pi / 2  # π/2


class ButterworthFilter(AnalogBasedFilter):
    """Базовый класс фильтра Баттерворта с вычислением нормированных полюсов."""

    def __init__(self, B: np.ndarray, A: np.ndarray, spec: Specification):
        super().__init__(B, A, spec)

    @staticmethod
    def get_norm_poles(N: int, EpsP: float, W0: float = 1.0) -> list:
        """
        Нормированные полюса фильтра Баттерворта.

        Полюса размещаются равномерно на окружности радиуса alpha в левой
        полуплоскости p-плоскости.

        Параметры
        ---------
        N    : порядок фильтра
        EpsP : неравномерность АЧХ в полосе пропускания (EpsP = sqrt(1/Gp² - 1))
        W0   : масштабирующий множитель (нормирующая частота)
        """
        if N <= 0:
            raise ValueError(f"Число полюсов должно быть больше 0, получено N={N}")

        r = N % 2           # нечётность порядка
        alpha = W0 * EpsP ** (-1.0 / N)  # радиус окружности полюсов
        th0 = _pi05 / N     # угловой шаг между полюсами

        poles = []
        if r != 0:
            poles.append(complex(-alpha, 0.0))  # центральный полюс нечётного порядка
        for i in range(r, N, 2):
            angle = th0 * (i - r + N + 1)
            z = alpha * np.exp(1j * angle)
            poles.append(z)
            poles.append(z.conjugate())
        return poles

    @staticmethod
    def get_norm_poles_gp(N: int, Gp: float, W0: float = 1.0) -> list:
        """
        Нормированные полюса Баттерворта с явным указанием Gp.

        Параметры
        ---------
        N  : порядок фильтра
        Gp : коэффициент передачи на границе пропускания (0 < Gp ≤ 1)
        W0 : масштабирующий множитель
        """
        if N <= 0:
            raise ValueError(f"Число полюсов должно быть больше 0, получено N={N}")

        r = N % 2
        # alpha = W0 * (1/Gp² - 1)^{-0.5/N}
        alpha = W0 * (1.0 / (Gp * Gp) - 1.0) ** (-0.5 / N)
        th0 = _pi05 / N

        poles = []
        if r != 0:
            poles.append(complex(-alpha, 0.0))
        for i in range(r, N, 2):
            angle = th0 * (i - r + N + 1)
            z = alpha * np.exp(1j * angle)
            poles.append(z)
            poles.append(z.conjugate())
        return poles
