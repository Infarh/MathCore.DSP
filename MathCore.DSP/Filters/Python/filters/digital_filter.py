# Цифровой фильтр — базовый класс (порт C# DigitalFilter)
import numpy as np
from scipy import signal as sp_signal


class DigitalFilter:
    """Базовый класс цифрового фильтра с методами билинейного преобразования."""

    # ── Билинейное пред-искажение частот ──────────────────────────────────────

    @staticmethod
    def to_digital_frequency(analog_freq: float, dt: float) -> float:
        """Аналоговая частота прототипа → пред-искажённая частота для билинейного преобразования."""
        return np.tan(np.pi * analog_freq * dt) / (np.pi * dt)

    @staticmethod
    def to_analog_frequency(digital_freq: float, dt: float) -> float:
        """Обратное преобразование: пред-искажённая частота → аналоговая."""
        return np.arctan(np.pi * digital_freq * dt) / (np.pi * dt)

    @staticmethod
    def to_digital_frequency_w(analog_w: float, dt: float) -> float:
        """Аналоговая циклическая частота → пред-искажённая циклическая частота."""
        return np.tan(analog_w * dt / 2) / (dt / 2)

    @staticmethod
    def to_analog_frequency_w(digital_w: float, dt: float) -> float:
        """Пред-искажённая циклическая частота → аналоговая циклическая частота."""
        return np.arctan(digital_w * dt / 2) / (dt / 2)

    # ── Билинейное z-преобразование полюсов/нулей ────────────────────────────

    @staticmethod
    def to_z(p: complex, dt: float) -> complex:
        """Полюс/нуль из p-плоскости в z-плоскость (билинейное преобразование)."""
        w = 2.0 / dt
        return (w + p) / (w - p)

    @staticmethod
    def to_z_array(p_list, dt: float, W0: float = 1.0) -> np.ndarray:
        """Массив полюсов/нулей из p-плоскости в z-плоскость с масштабированием W0."""
        return np.array([DigitalFilter.to_z(p * W0, dt) for p in p_list])

    @staticmethod
    def from_z(z: complex, dt: float) -> complex:
        """Полюс/нуль из z-плоскости обратно в p-плоскость."""
        w = 2.0 / dt
        return w * (z - 1.0) / (z + 1.0)

    # ── Частотная характеристика ──────────────────────────────────────────────

    @staticmethod
    def frequency_response(A: np.ndarray, B: np.ndarray, freqs, fd: float = 1.0) -> np.ndarray:
        """АЧХ/ФЧХ фильтра H(f) для заданных частот (Гц), нормированных на fd."""
        w = 2.0 * np.pi * np.asarray(freqs, dtype=float) / fd
        _, h = sp_signal.freqz(B, A, worN=w)
        return h
