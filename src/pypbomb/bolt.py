"""
Functions for handling detonation tube bolt calculations.
"""

import warnings
from dataclasses import dataclass
from enum import Enum
from functools import cached_property

import numpy as np

from pypbomb import _validate
from pypbomb.data import _units, bolt
from pypbomb.data.bolt import ExternalThreadAnsi, InternalThreadAnsi


@dataclass(frozen=True)
class SafetyFactors:
    bolt: float
    "Bolt safety factor"
    plate: float
    "Safety factor for plate threads (pullout)"


def safety_factors(
    internal_pressure: float,
    window_area: float,
    num_bolts: int,
    thread_size: str,
    partial_thread_class: str,
    bolt_max_tensile: float,
    plate_max_tensile: float,
    engagement_length: float,
    external_pressure: float = 101_325.0,
) -> SafetyFactors:
    """
    Calculate bolt and plate safety factors for viewing window bolts

    :param internal_pressure: Pressure inside the tube (Pa)
    :param external_pressure: Pressure outside the tube (Pa) - defaults to 1 atm
    :param window_area: Exposed (unclamped) window area (m^2)
    :param num_bolts: Number of bolts attaching viewing window to tube
    :param thread_size: Size of bolt thread, e.g. "1/4-20" (ANSI inch)
    :param partial_thread_class: "2" for 2A/2B threads, "3" for 3A/3B
    :param bolt_max_tensile: Bolt (external thread) tensile failure stress (Pa)
    :param plate_max_tensile: Plate (internal thread) tensile failure stress (Pa)
    :param engagement_length: Total thread engagement length (m)
    :return: Safety factor for bolt and plate for the given pressure differential
    """
    areas = stress_areas(
        thread_size=thread_size,
        partial_thread_class=partial_thread_class,
        bolt_max_tensile=bolt_max_tensile,
        engagement_length=engagement_length,
    )
    return _safety_factors_from_stress_areas(
        areas=areas,
        internal_pressure=internal_pressure,
        window_area=window_area,
        num_bolts=num_bolts,
        bolt_max_tensile=bolt_max_tensile,
        plate_max_tensile=plate_max_tensile,
        external_pressure=external_pressure,
    )


class BoltFailureMode(Enum):
    tension = "tension"
    shear = "shear"


@dataclass(frozen=True)
class BoltStressAreas:
    shear: float
    "External thread (bolt) shear stress area (m^2)"
    tensile: float
    "External thread (bolt) tensile stress area (m^2)"
    failure_mode: BoltFailureMode
    "Limiting factor for bolt stress. Plate may be damaged in an overpressure event if bolts are shear-limited"

    @cached_property
    def limiting(self) -> float:
        """
        Limiting (minimum) stress area
        """
        return min(self.shear, self.tensile)


@dataclass(frozen=True)
class StressAreas:
    bolt: BoltStressAreas
    plate: float
    "Internal thread (plate) shear area (m^2)"

    def safety_factors(
        self,
        internal_pressure: float,
        window_area: float,
        num_bolts: int,
        bolt_max_tensile: float,
        plate_max_tensile: float,
        external_pressure: float = 101_325.0,
    ) -> SafetyFactors:
        """
        Calculate bolt and plate safety factors for viewing window bolts

        :param internal_pressure: Pressure inside the tube (Pa)
        :param external_pressure: Pressure outside the tube (Pa) - defaults to 1 atm
        :param window_area: Exposed (unclamped) window area (m^2)
        :param num_bolts: Number of bolts attaching viewing window to tube
        :param bolt_max_tensile: Bolt (external thread) tensile failure stress (Pa)
        :param plate_max_tensile: Plate (internal thread) tensile failure stress (Pa)
        :return: Safety factor for bolt and plate for the given pressure differential
        """
        return _safety_factors_from_stress_areas(
            areas=self,
            internal_pressure=internal_pressure,
            window_area=window_area,
            num_bolts=num_bolts,
            bolt_max_tensile=bolt_max_tensile,
            plate_max_tensile=plate_max_tensile,
            external_pressure=external_pressure,
        )


def _safety_factors_from_stress_areas(
    areas: StressAreas,
    internal_pressure: float,
    window_area: float,
    num_bolts: int,
    bolt_max_tensile: float,
    plate_max_tensile: float,
    external_pressure: float = 101_325.0,
):
    _validate.pressure(internal_pressure)
    _validate.pressure(external_pressure)
    _validate.area(window_area)
    _validate.pressure(bolt_max_tensile)
    _validate.pressure(plate_max_tensile)
    if num_bolts < 1:
        raise ValueError("At least one bolt is required")

    force_per_bolt = (internal_pressure - external_pressure) * window_area / num_bolts

    return SafetyFactors(
        bolt=bolt_max_tensile / (force_per_bolt / areas.bolt.limiting),
        plate=plate_max_tensile / (force_per_bolt / areas.plate),
    )


def stress_areas(
    thread_size: str,
    partial_thread_class: str,
    bolt_max_tensile: float,
    engagement_length: float,
) -> StressAreas:
    """
    Calculates internal and external thread stress areas using formulas from Machinery's Handbook, 26th edition.

    :param thread_size: Size of bolt thread, e.g. "1/4-20" (ANSI inch)
    :param partial_thread_class: "2" for 2A/2B threads, "3" for 3A/3B
    :param bolt_max_tensile: Bolt (external thread) tensile failure stress (Pa)
    :param engagement_length: Total thread engagement length (m)
    :return: c, along with predicted bolt failure mode
    """
    _validate.pressure(bolt_max_tensile)
    _validate.length(engagement_length)

    bolt_dimensions = bolt.external_thread_dimensions(thread_class=f"{partial_thread_class}A", thread_size=thread_size)
    plate_dimensions = bolt.internal_thread_dimensions(thread_class=f"{partial_thread_class}B", thread_size=thread_size)

    tpm = bolt_dimensions.threads_per_meter
    k_n_max = plate_dimensions.diameter.minor.max
    e_s_min = bolt_dimensions.diameter.pitch.min
    e_n_max = plate_dimensions.diameter.pitch.max
    d_s_min = bolt_dimensions.diameter.major.min
    d_basic = bolt_dimensions.diameter.basic

    bolt_area_shear = _bolt_area_shear(tpm=tpm, engagement_length=engagement_length, k_n_max=k_n_max, e_s_min=e_s_min)
    if bolt_max_tensile < _units.kpsi_to_pa(100):
        bolt_area_tensile = _bolt_area_tensile_lt_100kpsi(tpm=tpm, d_basic=d_basic)
    else:
        bolt_area_tensile = _bolt_area_tensile_ge_100kpsi(tpm=tpm, e_s_min=e_s_min)

    if bolt_area_shear < bolt_area_tensile:
        warnings.warn(
            "Bolts fail in shear, not tension. Plate may be damaged. Consider increasing bolt engagement length.",
            RuntimeWarning,
        )
        bolt_failure_mode = BoltFailureMode.shear
    else:
        bolt_failure_mode = BoltFailureMode.tension

    plate_area = _plate_area_shear(tpm=tpm, engagement_length=engagement_length, d_s_min=d_s_min, e_n_max=e_n_max)

    return StressAreas(
        bolt=BoltStressAreas(
            shear=bolt_area_shear,
            tensile=bolt_area_tensile,
            failure_mode=bolt_failure_mode,
        ),
        plate=plate_area,
    )


def _bolt_area_shear(tpm: float, engagement_length: float, k_n_max: float, e_s_min: float) -> float:
    """
    Calculate screw shear area using eq. 5 (p. 1491) in Fasteners section of Machinery's Handbook 26

    :param tpm: Threads per meter (1/m)
    :param engagement_length: Total thread engagement length (m)
    :param k_n_max: Maximum plate minor diameter (m)
    :param e_s_min: Minimum bolt pitch diameter (m)
    :return: Bolt shear area (m^2)
    """
    return np.pi * tpm * engagement_length * k_n_max * (1 / (2 * tpm) + 0.57735 * (e_s_min - k_n_max))


def _bolt_area_tensile_lt_100kpsi(tpm: float, d_basic: float) -> float:
    """
    Calculate screw tensile area for tensile strengths < 100 kpsi using eqs. 9 (p. 1482) and 2a (p. 1490) in Fasteners
    section of Machinery's Handbook 26

    :param tpm: Threads per meter (1/m)
    :param d_basic: Bolt basic diameter (m)
    :return: Bolt tensile area (m^2)
    """
    return np.pi / 4 * (d_basic - 0.9742785 / tpm) ** 2


def _bolt_area_tensile_ge_100kpsi(tpm: float, e_s_min: float) -> float:
    """
    Calculate screw tensile area for tensile strengths >= 100 kpsi using eq. 2b (p. 1490) in Fasteners section of
    Machinery's Handbook 26

    :param tpm: Threads per meter (1/m)
    :param e_s_min: Minimum bolt pitch diameter (m)
    :return: Bolt tensile area (m^2)
    """
    return np.pi * (e_s_min / 2 - 0.16238 / tpm) ** 2


def _plate_area_shear(tpm: float, engagement_length: float, d_s_min: float, e_n_max: float) -> float:
    """
    Calculate plate shear area using eq. 6 (p. 1491) in Fasteners section of Machinery's Handbook 26

    :param tpm: Threads per meter (1/m)
    :param engagement_length: Total thread engagement length (m)
    :param d_s_min: Minimum bolt major diameter (m)
    :param e_n_max: Maximum plate pitch diameter
    :return: Plate shear area (m^2)
    """
    return np.pi * tpm * engagement_length * d_s_min * (1 / (2 * tpm) + 0.57735 * (d_s_min - e_n_max))


def minimum_engagement(
    bolt_dimensions: ExternalThreadAnsi,
    plate_dimensions: InternalThreadAnsi,
    areas: StressAreas,
    bolt_max_tensile: float,
    plate_max_tensile: float,
) -> float:
    """
    Calculate minimum thread engagement (corrected for material differences) using eqs. 1 and 4 (pp. 1490-1491) in
    Fasteners section of Machinery's Handbook 26

    :param bolt_dimensions: Bolt dimensions (ANSI external thread)
    :param plate_dimensions: Plate dimensions (ANSI internal thread)
    :param areas: Bolt and plate stress areas
    :param bolt_max_tensile: Bolt max tensile strength (Pa)
    :param plate_max_tensile: Plate max tensile strength (Pa)
    :return: Minimum engagement length causing bolt to fail in tension rather than shear, thus preventing the plate
        from stripping during an overpressure event
    """
    _validate.pressure(bolt_max_tensile)
    _validate.pressure(plate_max_tensile)

    tpm = bolt_dimensions.threads_per_meter
    k_n_max = plate_dimensions.diameter.minor.max
    e_s_min = bolt_dimensions.diameter.pitch.min

    j = _j_factor(
        areas=areas,
        bolt_max_tensile=bolt_max_tensile,
        plate_max_tensile=plate_max_tensile,
    )

    return (2 * areas.bolt.tensile / (k_n_max * np.pi * (0.5 + 0.57735 * tpm * (e_s_min - k_n_max)))) * j


def _j_factor(areas: StressAreas, bolt_max_tensile: float, plate_max_tensile: float) -> float:
    """
    Calculate J-factor using eq.3 (p. 1490) from Fasteners section of Machinery's Handbook 26

    :param areas: Bolt and plate stress areas
    :param bolt_max_tensile: Bolt max tensile strength (Pa)
    :param plate_max_tensile: Plate max tensile strength (Pa)
    :return: J-factor
    """
    return (areas.bolt.shear * bolt_max_tensile) / (areas.plate * plate_max_tensile)
