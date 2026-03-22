# Базовый класс эллиптического фильтра (порт C# EllipticFilter)
# Включает реализацию специальных функций Якоби с комплексными аргументами
import numpy as np
from math import ceil

try:
    import mpmath
    _HAS_MPMATH = True
except ImportError:
    _HAS_MPMATH = False

from scipy.special import ellipk, ellipj

from .analog_based_filter import AnalogBasedFilter, Specification


# ── Специальные эллиптические функции ─────────────────────────────────────────

def _ellipk(k: float) -> float:
    """Полный эллиптический интеграл первого рода K(k). Аргумент — модуль k, не k²."""
    return float(ellipk(k * k))


def _ellipk_complementary(k: float) -> float:
    """Дополнительный эллиптический интеграл K'(k) = K(√(1-k²))."""
    return float(ellipk(1.0 - k * k))


def sn_uk_real(u: float, k: float) -> float:
    """Нормированная функция Якоби sn(u·K(k), k) для вещественного u."""
    m = k * k
    K = float(ellipk(m))
    sn, cn, dn, _ = ellipj(u * K, m)
    return float(sn)


def cd_uk_real(u: float, k: float) -> float:
    """Нормированная функция Якоби cd(u·K(k), k) = cn/dn для вещественного u."""
    m = k * k
    K = float(ellipk(m))
    sn, cn, dn, _ = ellipj(u * K, m)
    return float(cn / dn)


def sn_uk(u, k: float):
    """
    Нормированная функция Якоби sn(u·K(k), k) для вещественного или комплексного u.

    Для комплексных аргументов использует mpmath. Для вещественных — scipy.
    """
    if isinstance(u, complex) and u.imag != 0.0:
        if not _HAS_MPMATH:
            raise ImportError("mpmath требуется для вычисления функций Якоби с комплексным аргументом. "
                              "Установите: pip install mpmath")
        m = float(k) ** 2
        K_val = mpmath.ellipk(m)
        u_mp = mpmath.mpc(u.real, u.imag) * K_val
        return complex(mpmath.ellipfun('sn', u_mp, m))
    return sn_uk_real(float(u.real if isinstance(u, complex) else u), k)


def cd_uk(u, k: float):
    """
    Нормированная функция Якоби cd(u·K(k), k) для вещественного или комплексного u.

    Для комплексных аргументов использует mpmath. Для вещественных — scipy.
    """
    if isinstance(u, complex) and u.imag != 0.0:
        if not _HAS_MPMATH:
            raise ImportError("mpmath требуется для вычисления функций Якоби с комплексным аргументом. "
                              "Установите: pip install mpmath")
        m = float(k) ** 2
        K_val = mpmath.ellipk(m)
        u_mp = mpmath.mpc(u.real, u.imag) * K_val
        cn_val = mpmath.ellipfun('cn', u_mp, m)
        dn_val = mpmath.ellipfun('dn', u_mp, m)
        return complex(cn_val / dn_val)
    return cd_uk_real(float(u.real if isinstance(u, complex) else u), k)


def sn_inverse(w: complex, k: float) -> complex:
    """
    Нормированный обратный синус Якоби sn⁻¹(w, k) / K(k) для комплексного w.

    Вычисляет F(arcsin(w), k²) / K(k) через mpmath.
    """
    if not _HAS_MPMATH:
        raise ImportError("mpmath требуется для sn_inverse. Установите: pip install mpmath")
    m = float(k) ** 2
    K_val = mpmath.ellipk(m)
    w_mp = mpmath.mpc(w.real, w.imag)
    phi = mpmath.asin(w_mp)
    F_val = mpmath.ellipf(phi, m)
    return complex(F_val / K_val)


# ── Базовый класс ──────────────────────────────────────────────────────────────

class EllipticFilter(AnalogBasedFilter):
    """
    Базовый класс эллиптического фильтра.

    Реализует вычисление нормированных нулей и полюсов эллиптического ФНЧ
    на основе функций Якоби и эллиптических интегралов.
    Имеет наилучшую крутизну АЧХ при заданных ограничениях на пульсации.
    """

    def __init__(self, B: np.ndarray, A: np.ndarray, spec: Specification):
        super().__init__(B, A, spec)

    # ── Эллиптические интегралы ───────────────────────────────────────────────

    @staticmethod
    def K(k: float) -> float:
        """Полный эллиптический интеграл первого рода K(k)."""
        return _ellipk(k)

    @staticmethod
    def T(k: float) -> float:
        """Дополнительный эллиптический интеграл K'(k) = K(√(1-k²))."""
        return _ellipk_complementary(k)

    # ── Нули эллиптического прототипа ─────────────────────────────────────────

    @staticmethod
    def enum_normed_zeros(N: int, EpsP: float, EpsS: float, W0: float = 1.0) -> list:
        """
        Нормированные нули эллиптического ФНЧ-прототипа.

        Нули расположены на мнимой оси в точках ±j·W0/(k_w·cd_uk(u_i, k_w)).

        Параметры
        ---------
        N    : порядок фильтра
        EpsP : неравномерность в полосе пропускания
        EpsS : неравномерность в полосе заграждения
        W0   : нормирующий множитель частоты
        """
        k_eps = EpsP / EpsS
        L = N // 2

        u = [(2 * i - 1.0) / N for i in range(1, L + 1)]

        m = np.sqrt(1.0 - k_eps ** 2)
        kp = m ** N * np.prod([sn_uk_real(ui, m) ** 4 for ui in u]) if L > 0 else m ** N
        k_w = np.sqrt(1.0 - kp ** 2)

        zeros = []
        for ui in u:
            cd_val = cd_uk_real(ui, k_w)
            im_val = W0 / (k_w * cd_val)
            zeros.append(complex(0.0, +im_val))
            zeros.append(complex(0.0, -im_val))
        return zeros

    # ── Полюса эллиптического прототипа ───────────────────────────────────────

    @staticmethod
    def enum_normed_poles(N: int, EpsP: float, EpsS: float, W0: float = 1.0) -> list:
        """
        Нормированные полюса эллиптического ФНЧ-прототипа.

        Вычисляются через функции Якоби cd и sn с комплексным аргументом v0.
        """
        k_eps = EpsP / EpsS
        L, r = divmod(N, 2)

        u = [(2 * i - 1.0) / N for i in range(1, L + 1)]

        m = np.sqrt(1.0 - k_eps ** 2)
        kp = m ** N * np.prod([sn_uk_real(ui, m) ** 4 for ui in u]) if L > 0 else m ** N
        k_w = np.sqrt(1.0 - kp ** 2)
        v0 = sn_inverse(complex(0.0, 1.0 / EpsP), k_eps) / N

        poles = []
        if r != 0:
            # нечётный полюс: j·W0·sn(v0, k_w)
            sn_val = sn_uk(v0, k_w)
            poles.append(complex(0.0, W0) * sn_val)

        for ui in u:
            # cd(u_i - v0, k_w): вещественная часть → p_im, мнимая → p_re
            cd_val = cd_uk(ui - v0, k_w) * W0
            p_im = cd_val.real   # соответствует (p_im, p_re) в C# Deconstruct
            p_re = cd_val.imag
            poles.append(complex(-p_re, +p_im))
            poles.append(complex(-p_re, -p_im))

        return poles

    # ── Нули и полюса вместе ──────────────────────────────────────────────────

    @staticmethod
    def get_normed_zeros_poles(N: int, EpsP: float, EpsS: float, W0: float = 1.0):
        """
        Нормированные нули и полюса эллиптического ФНЧ-прототипа.

        Возвращает (zeros, poles) — два списка комплексных чисел.
        Нулей всегда чётное число (N - r), полюсов — N.
        """
        k_eps = EpsP / EpsS
        L, r = divmod(N, 2)

        u = [(2 * i - 1.0) / N for i in range(1, L + 1)]

        m = np.sqrt(1.0 - k_eps ** 2)
        kp = m ** N * np.prod([sn_uk_real(ui, m) ** 4 for ui in u]) if L > 0 else m ** N
        k_w = np.sqrt(1.0 - kp ** 2)
        v0 = sn_inverse(complex(0.0, 1.0 / EpsP), k_eps) / N

        poles = [None] * N
        zeros = [None] * (N - r)

        if r != 0:
            sn_val = sn_uk(v0, k_w)
            poles[0] = complex(0.0, W0) * sn_val

        for i, ui in enumerate(u):
            cd_val = cd_uk(ui - v0, k_w) * W0
            p_im = cd_val.real   # C# deconstruct: (p_im, p_re) = cd_val
            p_re = cd_val.imag
            poles[2 * i + r]     = complex(-p_re, +p_im)
            poles[2 * i + r + 1] = complex(-p_re, -p_im)

            im_z = W0 / (k_w * cd_uk_real(ui, k_w))
            zeros[2 * i]     = complex(0.0, +im_z)
            zeros[2 * i + 1] = complex(0.0, -im_z)

        return zeros, poles
