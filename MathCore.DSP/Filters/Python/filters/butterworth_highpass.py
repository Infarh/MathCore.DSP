# Высокочастотный фильтр Баттерворта (порт C# ButterworthHighPass)
import numpy as np
from math import ceil, log

from .butterworth_filter import ButterworthFilter
from .analog_based_filter import AnalogBasedFilter, Specification
from .digital_filter import DigitalFilter


class ButterworthHighPass(ButterworthFilter):
    """
    Цифровой высокочастотный фильтр Баттерворта (ФВЧ).

    Проектируется через ФНЧ-прототип Баттерворта с LP→HP трансформацией
    и последующим билинейным z-преобразованием.
    """

    @staticmethod
    def get_polynoms(spec: Specification):
        """
        Коэффициенты передаточной функции ФВЧ Баттерворта.

        Алгоритм
        --------
        1. N = ceil(log(kEps) / (-log(kW))).
        2. Нормированные полюса Баттерворта.
        3. LP→HP трансформация: p → Wp/p.
        4. Билинейное z-преобразование.
        5. Нормирующий коэффициент g = |prod((1+zi)/2)| (Найквист z=-1).
        6. Числитель B — знакочередующиеся биномиальные коэффициенты.
        7. Знаменатель A = np.poly(z_poles).real.
        """
        N = int(ceil(log(spec.kEps) / (-log(spec.kW))))

        poles = ButterworthFilter.get_norm_poles(N, spec.EpsP)

        # LP→HP трансформация
        hp_poles = AnalogBasedFilter.transform_to_high_pass_w(poles, spec.Wp)

        # билинейное z-преобразование
        z_poles = DigitalFilter.to_z_array(hp_poles, spec.dt)

        # нормирующий коэффициент по Найквисту (z = -1)
        g_norm = abs(np.prod([(1.0 + z) / 2.0 for z in z_poles]))

        # числитель: знакочередующиеся биномиальные коэффициенты
        binomial = AnalogBasedFilter.get_binomial_coefficients(N)
        B = np.array([b * (g_norm if i % 2 == 0 else -g_norm)
                      for i, b in enumerate(binomial)])

        A = np.poly(z_poles).real

        return A, B

    def __init__(self, dt: float, fs: float, fp: float,
                 Gp: float = 0.891250938, Gs: float = 0.01):
        """
        Создаёт ФВЧ Баттерворта.

        Параметры
        ---------
        dt : период дискретизации (с)
        fs : частота среза полосы заграждения (Гц), fs < fp
        fp : граничная частота полосы пропускания (Гц)
        Gp : коэффициент передачи на границе пропускания
        Gs : коэффициент передачи на границе заграждения
        """
        spec = AnalogBasedFilter.get_specification(dt, fs, fp, Gp, Gs)
        A, B = ButterworthHighPass.get_polynoms(spec)
        super().__init__(B, A, spec)
