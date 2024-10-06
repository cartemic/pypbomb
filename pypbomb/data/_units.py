from typing import Optional

import numpy as np
from numpy import nan

from pypbomb._types import FloatArray

_convertable = float | FloatArray
_convertable_opt = Optional[_convertable]


def _scale_scalar_gt0(n: float, factor: float) -> float:
    return nan if n < 0 else n * factor


def _scale_array_gt0(n: FloatArray, factor: float) -> FloatArray:
    # FloatArray doesn't seem to warn about int values, which raise hell when being replaced by nan
    scaled = n.astype(float) * factor
    scaled[scaled < 0] = nan
    return scaled


def _convert(n: _convertable_opt, factor: float):
    if n is None:
        return nan
    if isinstance(n, np.ndarray):
        return _scale_array_gt0(n, factor)
    return _scale_scalar_gt0(n, factor)


def in_to_m(n: _convertable_opt) -> _convertable:
    return _convert(n, 0.0254)


def tpi_to_tpm(n: _convertable_opt) -> _convertable:
    return _convert(n, 0.0254**-1)


def gpa_to_pa(n: _convertable_opt) -> _convertable:
    return _convert(n, 1_000_000_000)


def kpsi_to_pa(n: _convertable_opt) -> _convertable:
    return _convert(n, 6_894_757.2931783)


def psi_to_pa(n: _convertable_opt) -> _convertable:
    return _convert(n, 6894.7572931783)


def gcm3_to_kgm3(n: _convertable_opt) -> _convertable:
    return _convert(n, 1_000)


def k_to_f(n: _convertable_opt) -> _convertable:
    if n is None:
        return nan
    elif isinstance(n, np.ndarray):
        scaled = ((n - 273.15) * 9 / 5) + 32
        scaled[n < 0] = nan
        return scaled
    elif n < 0:
        return nan
    return ((n - 273.15) * 9 / 5) + 32


def f_to_k(n: _convertable_opt) -> _convertable:
    if n is None:
        return nan
    elif isinstance(n, np.ndarray):
        scaled = ((n - 32) * 5 / 9) + 273.15
        scaled[scaled < 0] = nan
        return scaled
    scaled = ((n - 32) * 5 / 9) + 273.15
    if scaled < 0:
        return nan
    return scaled
