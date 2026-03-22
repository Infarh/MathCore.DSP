# БИХ-фильтр (IIR) — базовый класс (порт C# IIR)
import numpy as np
from scipy import signal as sp_signal
from dataclasses import dataclass
from typing import Optional, List, Tuple

from .digital_filter import DigitalFilter

_SOS_ORDER_THRESHOLD = 32


@dataclass
class _SosSection:
    """Секция второго порядка для SOS-оптимизации БИХ-фильтра."""
    b0: float
    b1: float
    b2: float
    a1: float
    a2: float


class IIR(DigitalFilter):
    """
    Фильтр с бесконечной импульсной характеристикой (БИХ).

    Реализует рекуррентное уравнение:
        A[0]*y[n] = B[0]*x[n] + B[1]*x[n-1] + ... - A[1]*y[n-1] - ...

    При порядке >= _SOS_ORDER_THRESHOLD автоматически строится
    SOS-разложение для улучшения численной стабильности.
    """

    def __init__(self, B: np.ndarray, A: np.ndarray):
        if len(B) == 0:
            raise ValueError("Массив коэффициентов числителя B не может быть пустым")
        if len(A) < 2:
            raise ValueError("Массив коэффициентов знаменателя A должен содержать не менее 2 элементов")
        if len(B) > len(A):
            raise ValueError("Число коэффициентов B не должно превышать число коэффициентов A")

        self._B = np.asarray(B, dtype=float)
        self._A = np.asarray(A, dtype=float)
        self._state: np.ndarray = np.zeros(max(len(A), len(B)) - 1)

        # SOS-оптимизация для высоких порядков
        self._sos_sections: Optional[List[_SosSection]] = None
        self._sos_gain: float = 1.0
        if self.order >= _SOS_ORDER_THRESHOLD:
            ok, sections, gain = IIR._try_create_sos(self._A, self._B)
            if ok:
                self._sos_sections = sections
                self._sos_gain = gain

    # ── Свойства ──────────────────────────────────────────────────────────────

    @property
    def B(self) -> np.ndarray:
        """Коэффициенты полинома числителя передаточной функции."""
        return self._B.copy()

    @property
    def A(self) -> np.ndarray:
        """Коэффициенты полинома знаменателя передаточной функции."""
        return self._A.copy()

    @property
    def order(self) -> int:
        """Порядок фильтра."""
        return max(len(self._A), len(self._B)) - 1

    # ── Обработка сигнала ─────────────────────────────────────────────────────

    def process(self, samples) -> np.ndarray:
        """Фильтрует массив отсчётов (без сохранения состояния)."""
        samples = np.asarray(samples, dtype=float)
        if self._sos_sections is not None:
            return sp_signal.sosfilt(self._to_scipy_sos(), samples)
        return sp_signal.lfilter(self._B, self._A, samples)

    def process_sample(self, sample: float) -> float:
        """Фильтрует один отсчёт с сохранением внутреннего состояния."""
        y, self._state = sp_signal.lfilter(self._B, self._A, [sample], zi=self._state)
        return float(y[0])

    def reset(self):
        """Сбрасывает внутреннее состояние фильтра."""
        self._state = np.zeros(max(len(self._A), len(self._B)) - 1)

    # ── АЧХ/ФЧХ ──────────────────────────────────────────────────────────────

    def get_frequency_response(self, freqs, fd: float = 1.0) -> np.ndarray:
        """АЧХ/ФЧХ H(f) для заданных частот в Гц."""
        return DigitalFilter.frequency_response(self._A, self._B, freqs, fd)

    # ── SOS-разложение ────────────────────────────────────────────────────────

    def get_sos(self) -> np.ndarray:
        """Возвращает SOS-представление в формате scipy [[b0,b1,b2,a0,a1,a2], ...]."""
        if self._sos_sections is not None:
            return self._to_scipy_sos()
        return sp_signal.tf2sos(self._B, self._A)

    def _to_scipy_sos(self) -> np.ndarray:
        """Конвертирует внутренние SOS-секции в формат scipy."""
        sos = []
        for s in self._sos_sections:
            sos.append([s.b0, s.b1, s.b2, 1.0, s.a1, s.a2])
        return np.array(sos)

    # ── Последовательное/параллельное соединение ──────────────────────────────

    def serial(self, other: 'IIR') -> 'IIR':
        """Последовательное соединение двух фильтров (свёртка числителей и знаменателей)."""
        return IIR(
            B=np.polymul(self._B, other._B),
            A=np.polymul(self._A, other._A),
        )

    def parallel(self, other: 'IIR') -> 'IIR':
        """Параллельное соединение двух фильтров."""
        return IIR(
            B=np.polymul(
                np.polyadd(np.polymul(self._B, other._A), np.polymul(other._B, self._A)),
                [0.5],
            ),
            A=np.polymul(self._A, other._A),
        )

    # ── Статические методы построения SOS ────────────────────────────────────

    @staticmethod
    def _try_create_sos(A: np.ndarray, B: np.ndarray) -> Tuple[bool, Optional[List[_SosSection]], float]:
        """Пытается построить SOS-разложение из коэффициентов передаточной функции."""
        try:
            a0 = A[0]
            if abs(a0) <= np.finfo(float).eps:
                return False, None, 1.0

            a_norm = A / a0
            b_norm = B / a0

            poles = np.roots(a_norm)
            zeros = np.roots(b_norm)

            den_sections = IIR._build_sections_from_roots(poles)
            if len(den_sections) == 0:
                return False, None, 1.0

            num_sections = IIR._build_sections_from_roots(zeros)
            if len(num_sections) > len(den_sections):
                return False, None, 1.0

            sections: List[_SosSection] = []
            for i, (_, a1, a2) in enumerate(den_sections):
                b0, b1, b2 = num_sections[i] if i < len(num_sections) else (1.0, 0.0, 0.0)
                sections.append(_SosSection(b0, b1, b2, a1, a2))

            gain = float(b_norm[0])
            IIR._distribute_gain(sections, gain)

            if not IIR._is_valid_sos(sections):
                return False, None, 1.0

            return True, sections, 1.0
        except Exception:
            return False, None, 1.0

    @staticmethod
    def _build_sections_from_roots(roots: np.ndarray) -> List[Tuple[float, float, float]]:
        """Строит список SOS-секций (b0, b1, b2) из корней полинома."""
        if len(roots) == 0:
            return []

        imag_eps = 1e-9
        complex_roots = [r for r in roots if abs(r.imag) > imag_eps]
        real_roots = sorted(
            [r.real for r in roots if abs(r.imag) <= imag_eps],
            key=lambda x: -abs(x),
        )

        sections: List[Tuple[float, float, float]] = []

        # Комплексно-сопряжённые пары → секции 2-го порядка
        while complex_roots:
            r1 = complex_roots.pop(0)
            conjugate = r1.conjugate()
            pair_idx, min_delta = -1, float('inf')
            for i, r in enumerate(complex_roots):
                d = abs(r - conjugate)
                if d < min_delta:
                    min_delta = d
                    pair_idx = i
            if pair_idx < 0:
                real_roots.append(r1.real)
                continue
            r2 = complex_roots.pop(pair_idx)
            s, p = (r1 + r2).real, (r1 * r2).real
            sections.append((1.0, -s, p))

        # Вещественные корни → секции (по паре)
        for i in range(0, len(real_roots), 2):
            r1 = real_roots[i]
            r2 = real_roots[i + 1] if i + 1 < len(real_roots) else 0.0
            sections.append((1.0, -(r1 + r2), r1 * r2))

        return sections

    @staticmethod
    def _distribute_gain(sections: List[_SosSection], gain: float):
        """Равномерно распределяет усиление gain по всем секциям."""
        if not sections or gain == 0 or not np.isfinite(gain):
            return
        section_gain = abs(gain) ** (1.0 / len(sections))
        if not np.isfinite(section_gain) or section_gain == 0:
            return
        sign = np.sign(gain)
        for i, s in enumerate(sections):
            scale = section_gain * (sign if i == 0 else 1.0)
            sections[i] = _SosSection(s.b0 * scale, s.b1 * scale, s.b2 * scale, s.a1, s.a2)

    @staticmethod
    def _is_valid_sos(sections: List[_SosSection], max_abs: float = 1e6, max_pole_radius: float = 1.0001) -> bool:
        """Проверяет численную корректность SOS-разложения."""
        for s in sections:
            vals = [s.b0, s.b1, s.b2, s.a1, s.a2]
            if any(not np.isfinite(v) for v in vals) or any(abs(v) > max_abs for v in vals):
                return False
            disc = s.a1 ** 2 - 4 * s.a2
            if disc >= 0:
                sq = np.sqrt(disc)
                if abs((-s.a1 + sq) / 2) > max_pole_radius or abs((-s.a1 - sq) / 2) > max_pole_radius:
                    return False
            else:
                re = -s.a1 / 2
                im = np.sqrt(-disc) / 2
                if np.sqrt(re ** 2 + im ** 2) > max_pole_radius:
                    return False
        return True

    def __repr__(self) -> str:
        return f"IIR(order={self.order})"
