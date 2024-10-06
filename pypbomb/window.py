"""
Functions to assist with the design of viewing windows for optical access.
"""

import warnings

import numpy as np

from pypbomb import _validate


def safety_factor(
    length: float,
    width: float,
    thickness: float,
    internal_pressure: float,
    rupture_modulus: float,
    external_pressure: float = 101_325.0,
) -> float:
    """
    Calculate the safety factor of a clamped rectangular window at a given pressure.

    :param length: Window unsupported (viewing) length (m)
    :param width: Window unsupported (viewing) width (m)
    :param thickness: Window thickness (m)
    :param internal_pressure: Pressure inside the tube (Pa)
    :param external_pressure: Pressure outside the tube (Pa) - defaults to 1 atm
    :param rupture_modulus: Rupture modulus of desired window material
    :return: Window factor of safety
    """
    _validate.length(length)
    _validate.length(width)
    _validate.length(thickness)
    _validate.pressure(internal_pressure)
    _validate.pressure(external_pressure)
    _validate.rupture_modulus(rupture_modulus)

    if np.isclose(internal_pressure, external_pressure):
        return np.inf

    # https://www.crystran.co.uk/userfiles/files/design-of-pressure-windows.pdf
    window_pressure = abs(internal_pressure - external_pressure)
    k = 0.75  # clamped
    r = length / width
    return (thickness / length) ** 2 * (2 * rupture_modulus * (1 + r**2)) / (k * window_pressure)


def minimum_thickness(
    length: float,
    width: float,
    fs: float,
    internal_pressure: float,
    rupture_modulus: float,
    external_pressure: float = 101_325.0,
) -> float:
    """
    Calculate the thickness of a clamped rectangular window which gives the desired safety factor

    :param length: Window unsupported (viewing) length (m)
    :param width: Window unsupported (viewing) width (m)
    :param fs: Desired factor of safety
    :param internal_pressure: Pressure inside the tube (Pa)
    :param external_pressure: Pressure outside the tube (Pa) - defaults to 1 atm
    :param rupture_modulus:
    :return: Minimum window thickness required to achieve at least the desired factor of safety
    """
    _validate.length(length)
    _validate.length(width)
    _validate.rupture_modulus(rupture_modulus)
    if fs < 0:
        raise ValueError("Safety factor must be positive")
    elif fs < 1:
        warnings.warn(f"Window safety factor <1: {fs}", UserWarning)

    # https://www.crystran.co.uk/userfiles/files/design-of-pressure-windows.pdf
    window_pressure = abs(internal_pressure - external_pressure)
    k = 0.75  # clamped
    r = length / width
    thickness = length * np.sqrt(fs * k / 2) * np.sqrt(window_pressure / (rupture_modulus * (1 + r**2)))

    return thickness
