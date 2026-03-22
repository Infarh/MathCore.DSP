# Пакет Python-портов IIR-фильтров из MathCore.DSP
"""
filters — Python-порт цифровых IIR-фильтров библиотеки MathCore.DSP.

Иерархия классов
----------------
DigitalFilter
└── IIR
    └── AnalogBasedFilter
        ├── ButterworthFilter
        │   ├── ButterworthLowPass
        │   ├── ButterworthHighPass
        │   ├── ButterworthBandPass
        │   └── ButterworthBandStop
        ├── ChebyshevFilter
        │   ├── ChebyshevLowPass
        │   ├── ChebyshevHighPass
        │   ├── ChebyshevBandPass
        │   └── ChebyshevBandStop
        └── EllipticFilter
            ├── EllipticLowPass
            ├── EllipticHighPass
            ├── EllipticBandPass
            └── EllipticBandStop
"""

from .digital_filter import DigitalFilter
from .iir import IIR
from .analog_based_filter import AnalogBasedFilter, Specification

from .butterworth_filter import ButterworthFilter
from .butterworth_lowpass import ButterworthLowPass
from .butterworth_highpass import ButterworthHighPass
from .butterworth_bandpass import ButterworthBandPass
from .butterworth_bandstop import ButterworthBandStop

from .chebyshev_filter import ChebyshevFilter, ChebyshevType
from .chebyshev_lowpass import ChebyshevLowPass
from .chebyshev_highpass import ChebyshevHighPass
from .chebyshev_bandpass import ChebyshevBandPass
from .chebyshev_bandstop import ChebyshevBandStop

from .elliptic_filter import (
    EllipticFilter,
    sn_uk, cd_uk, sn_inverse,
    sn_uk_real, cd_uk_real,
)
from .elliptic_lowpass import EllipticLowPass
from .elliptic_highpass import EllipticHighPass
from .elliptic_bandpass import EllipticBandPass
from .elliptic_bandstop import EllipticBandStop

__all__ = [
    # базовые классы
    "DigitalFilter", "IIR", "AnalogBasedFilter", "Specification",
    # Баттерворт
    "ButterworthFilter",
    "ButterworthLowPass", "ButterworthHighPass",
    "ButterworthBandPass", "ButterworthBandStop",
    # Чебышев
    "ChebyshevFilter", "ChebyshevType",
    "ChebyshevLowPass", "ChebyshevHighPass",
    "ChebyshevBandPass", "ChebyshevBandStop",
    # Эллиптический
    "EllipticFilter",
    "sn_uk", "cd_uk", "sn_inverse", "sn_uk_real", "cd_uk_real",
    "EllipticLowPass", "EllipticHighPass",
    "EllipticBandPass", "EllipticBandStop",
]
