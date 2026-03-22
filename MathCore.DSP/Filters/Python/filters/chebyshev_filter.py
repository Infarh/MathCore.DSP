# Базовый класс фильтра Чебышева (порт C# ChebyshevFilter + ChebyshevType)
from enum import Enum

import numpy as np

from .analog_based_filter import AnalogBasedFilter, Specification


class ChebyshevType(Enum):
    """Тип (род) фильтра Чебышева."""
    I = 1             # равноволновая пульсация в полосе пропускания
    II = 2            # равноволновая пульсация в полосе заграждения
    IICorrected = 3   # II рода с точным попаданием в частоты среза и заграждения


class ChebyshevFilter(AnalogBasedFilter):
    """
    Базовый класс фильтра Чебышева.

    Содержит вычисление нормированных нулей и полюсов для обоих родов.
    """

    def __init__(self, B: np.ndarray, A: np.ndarray, spec: Specification,
                 filter_type: ChebyshevType):
        super().__init__(B, A, spec)
        self.filter_type = filter_type

    # ── Гиперболические обратные функции ──────────────────────────────────────

    @staticmethod
    def arcsh(x: float) -> float:
        """Арксинус гиперболический arcsh(x) = ln(x + √(x²+1))."""
        return np.log(x + np.sqrt(x * x + 1.0))

    @staticmethod
    def arcch(x: float) -> float:
        """Аркосинус гиперболический arcch(x) = ln(x + √(x²-1))."""
        return np.log(x + np.sqrt(x * x - 1.0))

    # ── Полюса Чебышева I рода ────────────────────────────────────────────────

    @staticmethod
    def get_normed_poles_i(N: int, EpsP: float, W0: float = 1.0) -> list:
        """
        Нормированные полюса ФНЧ Чебышева I рода.

        Полюса размещаются на эллипсе с полуосями sh (горизонтальная)
        и ch (вертикальная).

        Параметры
        ---------
        N    : порядок фильтра
        EpsP : неравномерность АЧХ в полосе пропускания
        W0   : нормирующий множитель частоты
        """
        beta = ChebyshevFilter.arcsh(1.0 / EpsP) / N
        sh = np.sinh(beta) * W0   # горизонтальная полуось эллипса
        ch = np.cosh(beta) * W0   # вертикальная полуось эллипса

        r = N % 2  # нечётность порядка
        poles = []
        if N % 2 != 0:
            poles.append(complex(-sh, 0.0))  # центральный вещественный полюс
        dth = np.pi / 2.0 / N
        for i in range(r, N, 2):
            th = dth * (i - r + 1)
            sin_th, cos_th = np.sin(th), np.cos(th)
            poles.append(complex(-sh * sin_th, +ch * cos_th))
            poles.append(complex(-sh * sin_th, -ch * cos_th))
        return poles

    # ── Полюса и нули Чебышева II рода ────────────────────────────────────────

    @staticmethod
    def get_normed_poles_ii(N: int, EpsS: float, W0: float = 1.0):
        """
        Нормированные нули и полюса ФНЧ Чебышева II рода.

        Нули расположены на мнимой оси (вне полосы пропускания),
        полюса — на деформированном эллипсе.

        Параметры
        ---------
        N    : порядок фильтра
        EpsS : неравномерность АЧХ в полосе заграждения
        W0   : нормирующий множитель частоты

        Возвращает
        ----------
        (zeros, poles) : кортеж комплексных массивов
        """
        beta = ChebyshevFilter.arcsh(EpsS) / N
        sh = np.sinh(beta)
        ch = np.cosh(beta)

        r = N % 2
        zeros = []
        poles = []

        if N % 2 != 0:
            poles.append(complex(-W0 / sh, 0.0))  # вещественный полюс нечётного порядка

        dth = np.pi / 2.0 / N
        for i in range(r, N, 2):
            th = dth * (i - r + 1)
            sin_th, cos_th = np.sin(th), np.cos(th)

            # промежуточный комплексный вектор
            z_re = -sh * sin_th
            z_im = ch * cos_th
            z_pow2 = z_re ** 2 + z_im ** 2  # |z|²
            norm = W0 / z_pow2

            # полюса: W0 / z.conjugate() = W0 * z.re / |z|² + j * W0 * (-z.im) / |z|²
            poles.append(complex(z_re * norm, +z_im * norm))
            poles.append(complex(z_re * norm, -z_im * norm))

            # нули на мнимой оси: ±j·W0/cos(th)
            zeros.append(complex(0.0, +W0 / cos_th))
            zeros.append(complex(0.0, -W0 / cos_th))

        return zeros, poles
