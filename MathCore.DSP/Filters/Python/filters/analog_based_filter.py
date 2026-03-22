# Базовый класс фильтра на основе аналогового прототипа (порт C# AnalogBasedFilter)
import numpy as np
from .iir import IIR
from .digital_filter import DigitalFilter

pi2 = 2.0 * np.pi  # 2π — часто используемая константа


class Specification:
    """
    Спецификация цифрового фильтра на основе аналогового прототипа.

    Хранит все параметры фильтра и вычисляет производные величины:
    частоты в аналоговой и цифровой областях, коэффициенты неравномерности АЧХ.
    """

    def __init__(self, dt: float, fp: float, fs: float, Gp: float, Gs: float):
        """
        Параметры
        ---------
        dt : период дискретизации (с)
        fp : граничная частота полосы пропускания (Гц)
        fs : граничная частота полосы заграждения (Гц)
        Gp : коэффициент передачи на границе пропускания (0 < Gs < Gp ≤ 1)
        Gs : коэффициент передачи на границе заграждения
        """
        if not (fp < 1.0 / (2.0 * dt)):
            raise ValueError("Частота пропускания fp должна быть меньше половины частоты дискретизации")
        if fp < 0:
            raise ValueError(f"Частота среза fp не может быть отрицательной: fp={fp}")
        if not (0 < Gs < Gp <= 1):
            raise ValueError(f"Требуется 0 < Gs < Gp ≤ 1, получено Gp={Gp}, Gs={Gs}")

        self.dt = dt
        self.fp = fp
        self.fs = fs
        self.Gp = Gp
        self.Gs = Gs

        # Затухания в дБ
        self.Rp = -20.0 * np.log10(Gp)   # потери в полосе пропускания
        self.Rs = -20.0 * np.log10(Gs)   # затухание в полосе заграждения

        # Неоднородности АЧХ (параметры аналогового прототипа)
        self.EpsP = np.sqrt(1.0 / (Gp * Gp) - 1.0)
        self.EpsS = np.sqrt(10.0 ** (self.Rs / 10.0) - 1.0)

        # Пред-искажённые частоты аналогового прототипа
        self.Fp = DigitalFilter.to_digital_frequency(fp, dt)
        self.Fs = DigitalFilter.to_digital_frequency(fs, dt)

        # Циклические частоты оригинального (аналогового) фильтра
        self.wp = pi2 * fp
        self.ws = pi2 * fs

        # Циклические частоты аналогового прототипа (пред-искажённые)
        self.Wp = pi2 * self.Fp
        self.Ws = pi2 * self.Fs

    @property
    def fd(self) -> float:
        """Частота дискретизации (Гц)."""
        return 1.0 / self.dt

    @property
    def kEps(self) -> float:
        """Отношение неоднородностей EpsS / EpsP."""
        return self.EpsS / self.EpsP

    @property
    def kw(self) -> float:
        """Отношение аналоговых частот fs / fp."""
        return self.fs / self.fp

    @property
    def kW(self) -> float:
        """Отношение цифровых частот Ws / Wp."""
        return self.Ws / self.Wp

    def __repr__(self) -> str:
        return (f"Specification(fd={self.fd:.1f} Гц, fp={self.fp:.1f} Гц, "
                f"fs={self.fs:.1f} Гц, Gp={self.Gp:.4f}, Gs={self.Gs:.4f})")


class AnalogBasedFilter(IIR):
    """
    Цифровой фильтр на основе аналогового прототипа.

    Реализует билинейное z-преобразование и трансформации нулей/полюсов
    нормированного ФНЧ в различные типы фильтров (ФНЧ, ФВЧ, ППФ, ПЗФ).
    """

    def __init__(self, B: np.ndarray, A: np.ndarray, spec: Specification):
        super().__init__(B, A)
        self.dt = spec.dt
        self.fp = spec.fp
        self.fs = spec.fs
        self.Gp = spec.Gp
        self.Gs = spec.Gs

    @property
    def fd(self) -> float:
        """Частота дискретизации (Гц)."""
        return 1.0 / self.dt

    @property
    def spec(self) -> Specification:
        """Спецификация фильтра."""
        return Specification(self.dt, self.fp, self.fs, self.Gp, self.Gs)

    # ── Статические методы получения спецификации ─────────────────────────────

    @staticmethod
    def get_specification(dt: float, fp: float, fs: float,
                          Gp: float = 0.891250938, Gs: float = 0.01) -> Specification:
        """Создаёт спецификацию фильтра."""
        return Specification(dt, fp, fs, Gp, Gs)

    # ── Трансформации нулей/полюсов ───────────────────────────────────────────

    @staticmethod
    def transform_to_low_pass_w(normed, wp: float) -> list:
        """ФНЧ-прототип → ФНЧ с частотой wp (в циклических частотах): p → p·wp."""
        return [p * wp for p in normed]

    @staticmethod
    def transform_to_high_pass_w(normed, wp: float) -> list:
        """ФНЧ-прототип → ФВЧ с частотой wp (в циклических частотах): p → wp/p."""
        return [wp / p for p in normed]

    @staticmethod
    def transform_to_band_pass_w(normed, w_min: float, w_max: float) -> list:
        """
        ФНЧ-прототип → ППФ (в циклических частотах).

        Каждый полюс/нуль p → пара (p·dW/2 ± sqrt((p·dW/2)² − Wc²)).
        """
        dw05 = (w_max - w_min) / 2.0
        wc2 = w_min * w_max
        result = []
        for p in normed:
            pdw = dw05 * p
            sqrt_val = np.sqrt(complex(pdw ** 2 - wc2))
            result.append(pdw + sqrt_val)
            result.append(pdw - sqrt_val)
        return result

    @staticmethod
    def transform_to_band_stop_w(normed, w_min: float, w_max: float) -> list:
        """
        ФНЧ-прототип → ПЗФ (в циклических частотах).

        Каждый полюс/нуль p → пара (dW/(2p) ± sqrt((dW/(2p))² − Wc²)).
        """
        dw05 = (w_max - w_min) / 2.0
        wc2 = w_min * w_max
        result = []
        for p in normed:
            pdw = dw05 / p
            sqrt_val = np.sqrt(complex(pdw ** 2 - wc2))
            result.append(pdw + sqrt_val)
            result.append(pdw - sqrt_val)
        return result

    @staticmethod
    def transform_to_low_pass(normed, fp: float) -> list:
        """ФНЧ-прототип → ФНЧ (в обычных частотах Гц)."""
        return AnalogBasedFilter.transform_to_low_pass_w(normed, pi2 * fp)

    @staticmethod
    def transform_to_high_pass(normed, fp: float) -> list:
        """ФНЧ-прототип → ФВЧ (в обычных частотах Гц)."""
        return AnalogBasedFilter.transform_to_high_pass_w(normed, pi2 * fp)

    @staticmethod
    def transform_to_band_pass(normed, fmin: float, fmax: float) -> list:
        """ФНЧ-прототип → ППФ (в обычных частотах Гц)."""
        return AnalogBasedFilter.transform_to_band_pass_w(normed, pi2 * fmin, pi2 * fmax)

    @staticmethod
    def transform_to_band_stop(normed, fmin: float, fmax: float) -> list:
        """ФНЧ-прототип → ПЗФ (в обычных частотах Гц)."""
        return AnalogBasedFilter.transform_to_band_stop_w(normed, pi2 * fmin, pi2 * fmax)

    # ── Вспомогательные методы ────────────────────────────────────────────────

    @staticmethod
    def get_binomial_coefficients(N: int) -> np.ndarray:
        """Биномиальные коэффициенты C(N,0), C(N,1), …, C(N,N) без целочисленного переполнения."""
        if N < 0:
            raise ValueError(f"Порядок бинома N не может быть отрицательным: N={N}")
        coeffs = np.zeros(N + 1, dtype=float)
        coeffs[0] = 1.0
        for k in range(1, N + 1):
            coeffs[k] = coeffs[k - 1] * (N - (k - 1)) / k
        return coeffs

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(fd={self.fd:.1f} Гц, fp={self.fp:.1f} Гц, fs={self.fs:.1f} Гц)"
