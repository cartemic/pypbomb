"""
Functions for designing and determining operational limits of a detonation tube.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Optional

from pypbomb import sd, thermochem
from pypbomb._types import SpeciesDefinition
from pypbomb.data.material import MaterialProperties
from pypbomb.data.piping import PipeDimensions


def dynamic_load_factor(
    pipe_dimensions: PipeDimensions,
    material_properties: MaterialProperties,
    cj_velocity: float,
    critical_velocity_band_frac: float = 0.1,
) -> float:
    """
    Calculate the dynamic load factor (DLF), by which a detonation tube's static analysis should be derated in order to
    account for the tube's response to transient pressures. DLF calculation is based on the work of Shepherd [1]. Since
    the limits of "approximately equal to" are not defined in the paper, a default value of plus or minus ten percent
    is assumed, thus ``plus_or_minus=0.1``.

    [1] Shepherd, J. E. (2009). *Structural Response of Piping to Internal Gas Detonation*. Journal of Pressure Vessel
    Technology, 131(3), 031204. https://doi.org/10.1115/1.3089497

    :param pipe_dimensions: Dimensions of piping used to build the tube
    :param material_properties: Material properties of piping used to build the tube
    :param cj_velocity: Chapman-Jouguet detonation velocity (m/s)
    :param critical_velocity_band_frac: Velocities falling within the tube's critical velocity +/- band % are considered
        approximately equal to the critical velocity for the purposes of dynamic load factor calculation. The default
        value of 0.1 means plus or minus ten percent. Must be between 0 and 1.
    :return: Dynamic load factor (unitless)
    """
    if not (0 < critical_velocity_band_frac < 1):
        raise ValueError("\nplus_or_minus factor not between 0 and 1")

    radius = pipe_dimensions.mean_diameter / 2
    crit_velocity = (
        (material_properties.elastic_modulus**2 * pipe_dimensions.wall_thickness**2)
        / (3 * material_properties.density**2 * radius**2 * (1 - material_properties.poisson**2))
    ) ** (1 / 4)

    # set limits for "approximately Vcrit"
    upper_bound = crit_velocity * (1.0 + critical_velocity_band_frac)
    lower_bound = crit_velocity * (1.0 - critical_velocity_band_frac)

    if cj_velocity < lower_bound:
        return 1.0
    elif cj_velocity > upper_bound:
        return 2.0
    else:
        return 4.0


def max_allowable_pressure(pipe_dimensions: PipeDimensions, max_stress: float, safety_factor: float = 1.0) -> float:
    """
    Calculate the maximum allowable pressure from the tube dimensions and stress limits using the basic longitudinal
    joint formula [1].

    [1] E. F. Megyesy, *Pressure vessel handbook*, Oklahoma City, OK: PV Publishing, Inc., 2001, p. 14.

    :param pipe_dimensions: Dimensions of piping used to build the tube
    :param max_stress: Maximum allowable stress in tube material, e.g. from ASME B31.1 (Pa)
    :param safety_factor: Desired safety factor
    :return: Pressure resulting in maximum allowable stress (Pa)
    """
    return max_stress * (2 * pipe_dimensions.wall_thickness) / (pipe_dimensions.mean_diameter * safety_factor)


@dataclass(frozen=True)
class MaxInitialPressureResult:
    pressure: float
    dlf: Optional[float]


class TubeKind(Enum):
    open = "open"
    closed = "closed"


def max_initial_pressure(
    max_pressure: float,
    initial_temperature: float,
    mechanism: str,
    species: SpeciesDefinition,
    pipe_dimensions: PipeDimensions,
    material_properties: MaterialProperties,
    initial_pressure_guess: float = 101325.0,
    dlf: Optional[float] = None,
    critical_velocity_band_frac: float = 0.1,
    rtol: float = 1e-4,
    max_iterations: int = 500,
    parallelize_cj_calc: bool = False,
    tube_kind: TubeKind = TubeKind.closed,
) -> MaxInitialPressureResult:
    """
    Calculate the reactant initial pressure which will result in the given maximum pressure within the specified
    detonation tube

    :param max_pressure: Maximum allowable pressure within the tube (Pa)
    :param initial_temperature: Mixture initial temperature (K)
    :param mechanism: Mechanism file name for Cantera
    :param species: Reactant species mole fractions
    :param pipe_dimensions: Dimensions of piping used to build the tube
    :param material_properties: Material properties of piping used to build the tube
    :param initial_pressure_guess: Starting guess for tube initial pressure (Pa)
    :param dlf: Prescribed dynamic load factor. If left as `None`, DLF will be calculated as needed.
    :param critical_velocity_band_frac: Velocities falling within the tube's critical velocity +/- band are considered
        approximately equal to the critical velocity for the purposes of dynamic load factor calculation. The default
        value of 0.1 means plus or minus ten percent. Must be between 0 and 1.
    :param rtol: Relative error tolerance below which initial pressure calculations are considered 'good enough'
    :param max_iterations: Maximum number of loop iterations before exit, defaults to 500
    :param parallelize_cj_calc: Use multiprocessing for CJ state calculation, which is faster but requires the function
        to be run from ``__main__``, which may not behave well when used via Jupyter
    :param tube_kind: Kind of tube (open or closed)
    :return: Initial pressure (Pa)
    """
    shock_pressure, cj_velocity = _shock_pressure_and_cj_velocity(
        initial_temperature=initial_temperature,
        initial_pressure=initial_pressure_guess,
        mechanism=mechanism,
        species=species,
        parallelize_cj_calc=parallelize_cj_calc,
        tube_kind=tube_kind,
    )
    if dlf is None:
        # calculate DLF as needed
        dlf_actual = dynamic_load_factor(
            pipe_dimensions=pipe_dimensions,
            material_properties=material_properties,
            cj_velocity=cj_velocity,
            critical_velocity_band_frac=critical_velocity_band_frac,
        )
    else:
        # use user-specified DLF
        dlf_actual = dlf

    iterations = 1  # initial calc above
    while _err(shock_pressure=shock_pressure, dlf=dlf_actual, max_pressure=max_pressure) > rtol:
        if iterations == max_iterations:
            raise StopIteration(f"Max initial pressure did not converge after {iterations} iterations")

        # update initial pressure guess
        initial_pressure_guess *= max_pressure / (dlf_actual * shock_pressure)

        # get reflected shock pressure
        shock_pressure, cj_velocity = _shock_pressure_and_cj_velocity(
            initial_temperature=initial_temperature,
            initial_pressure=initial_pressure_guess,
            mechanism=mechanism,
            species=species,
            parallelize_cj_calc=parallelize_cj_calc,
            tube_kind=tube_kind,
        )
        if dlf is None:
            dlf_actual = dynamic_load_factor(
                pipe_dimensions=pipe_dimensions,
                material_properties=material_properties,
                cj_velocity=cj_velocity,
                critical_velocity_band_frac=critical_velocity_band_frac,
            )

        iterations += 1

    return MaxInitialPressureResult(pressure=initial_pressure_guess, dlf=dlf_actual)


def _err(shock_pressure: float, dlf: float, max_pressure: float) -> float:
    return abs((shock_pressure * dlf - max_pressure) / max_pressure)


def _shock_pressure_and_cj_velocity(
    initial_temperature: float,
    initial_pressure: float,
    mechanism: str,
    species: SpeciesDefinition,
    parallelize_cj_calc: bool,
    tube_kind: TubeKind,
) -> tuple[float, float]:
    match tube_kind:
        case TubeKind.closed:
            shock = thermochem.reflected_cj_shock(
                initial_temperature=initial_temperature,
                initial_pressure=initial_pressure,
                species=species,
                mechanism=mechanism,
                parallelize_cj_calc=parallelize_cj_calc,
            )
            return shock.reflected.state.P, shock.cj.speed
        case TubeKind.open:
            shock = sd.cj.shock(
                initial_pressure=initial_pressure,
                initial_temperature=initial_temperature,
                species=species,
                mechanism=mechanism,
                parallelize=parallelize_cj_calc,
            )
            return shock.state.P, shock.speed
