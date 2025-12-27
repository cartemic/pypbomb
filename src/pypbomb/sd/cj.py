import multiprocessing as mp
import typing
import warnings
from dataclasses import dataclass

import cantera as ct
import numpy as np

from pypbomb._types import FloatArray, SpeciesDefinition
from pypbomb.sd import error, properties


@dataclass(frozen=True)
class CurveFit:
    a: float
    b: float
    c: float
    r2: float

    @staticmethod
    def empty() -> "CurveFit":
        return CurveFit(a=np.nan, b=np.nan, c=np.nan, r2=np.nan)


@dataclass(frozen=True)
class CjShock:
    speed: float
    speed_r2: float
    state: ct.Solution


def curve_fit(x: FloatArray, y: FloatArray) -> CurveFit:
    """
    Calculate the least squares fit of parabolic data. This is a vectorized version of
    ``sdtoolbox.PostShock.LSQ_CJspeed`` using ``np.linalg.lstsq``.
    """
    x = np.array((np.ones(len(x)), x, x**2)).T
    fit = np.linalg.lstsq(x, y, rcond=-1)
    thetas = fit[0]
    ss_res = fit[1][0]
    ss_tot = np.sum(np.square(y - np.mean(y)))
    # As of 2.1, numpy does not support shape typing, which hoses up type hinting. They recommend using `typing.cast`
    # to signal to the type checker that we know what we're doing [1,2].
    # [1] https://numpy.org/doc/2.1/reference/typing.html
    # [2] https://docs.python.org/3/library/typing.html#typing.cast
    return CurveFit(
        a=typing.cast(float, thetas[2]),
        b=typing.cast(float, thetas[1]),
        c=typing.cast(float, thetas[0]),
        r2=1 - ss_res / ss_tot,
    )


def _calc_cj_guess_at_density_ratio_and_set_state(
    working_gas: ct.Solution,
    initial_state_gas: ct.Solution,
    error_tol_temperature: float,
    error_tol_velocity: float,
    density_ratio: float,
    max_iterations: typing.Optional[int] = 500,
) -> float:
    """
    Set ``working_gas`` to the Chapman-Jouguet state, calculated using Reynolds' iterative method and return the

    This function corresponds to ``sdtoolbox.PostShock.CJ_calc``

    :param working_gas: Working gas mixture used for calculations.
    :param initial_state_gas: Working gas mixture in its initial, undetonated state.
    :param error_tol_temperature: Temperature error tolerance for iteration, in Kelvin.
    :param error_tol_velocity: Velocity error tolerance for iteration, in m/s.
    :param density_ratio: Density ratio.
    :param max_iterations: Maximum number of loop iterations used to calculate output.
    :return: Initial velocity resulting in the input density ratio, in m/s
    """
    initial_volume = 1 / initial_state_gas.density

    guess_temperature = 2000
    guess_velocity = 2000
    guess_volume = initial_volume / density_ratio
    guess_density = 1 / guess_volume
    properties.set_equilibrium(working_gas=working_gas, density=guess_density, temperature=guess_temperature)

    delta_temperature = 1000
    delta_velocity = 1000

    loop_counter = 0
    while abs(delta_temperature) > (error_tol_temperature * guess_temperature) or abs(delta_velocity) > (
        error_tol_velocity * guess_velocity
    ):
        loop_counter += 1
        # check for non-convergence
        if loop_counter == max_iterations:
            warnings.warn(
                "No convergence within {0} iterations".format(max_iterations),
                Warning,
            )
            return guess_velocity

        # calculate unperturbed enthalpy and press. error for current guess
        err_unperturbed = error.equilibrium(
            working_gas=working_gas,
            initial_state_gas=initial_state_gas,
            speed=guess_velocity,
        )

        # perturb temperature
        delta_temperature = 0.02 * guess_temperature
        perturbed_temperature = guess_temperature + delta_temperature
        properties.set_equilibrium(
            working_gas=working_gas,
            density=guess_density,
            temperature=perturbed_temperature,
        )

        # calculate error rates for temperature perturbed state
        err_perturbed = error.equilibrium(
            working_gas=working_gas,
            initial_state_gas=initial_state_gas,
            speed=guess_velocity,
        )
        derivative_enthalpy_temperature = (err_perturbed.enthalpy - err_unperturbed.enthalpy) / delta_temperature
        derivative_pressure_temperature = (err_perturbed.pressure - err_unperturbed.pressure) / delta_temperature

        # perturb velocity
        delta_velocity = 0.02 * guess_velocity
        perturbed_velocity = guess_velocity + delta_velocity
        perturbed_temperature = guess_temperature
        properties.set_equilibrium(
            working_gas=working_gas,
            density=guess_density,
            temperature=perturbed_temperature,
        )

        # calculate error rates for velocity perturbed state
        err_perturbed = error.equilibrium(
            working_gas=working_gas,
            initial_state_gas=initial_state_gas,
            speed=perturbed_velocity,
        )
        derivative_enthalpy_velocity = (err_perturbed.enthalpy - err_unperturbed.enthalpy) / delta_velocity
        derivative_pressure_velocity = (err_perturbed.pressure - err_unperturbed.pressure) / delta_velocity

        # invert matrix
        j = (
            derivative_enthalpy_temperature * derivative_pressure_velocity
            - derivative_pressure_temperature * derivative_enthalpy_velocity
        )
        b = [
            derivative_pressure_velocity,
            -derivative_enthalpy_velocity,
            -derivative_pressure_temperature,
            derivative_enthalpy_temperature,
        ]
        a = [-err_unperturbed.enthalpy, -err_unperturbed.pressure]

        delta_temperature = (b[0] * a[0] + b[1] * a[1]) / j
        delta_velocity = (b[2] * a[0] + b[3] * a[1]) / j

        # limit temperature changes
        max_temperature_delta = 0.2 * guess_temperature
        if abs(delta_temperature) > max_temperature_delta:
            delta_temperature *= max_temperature_delta / abs(delta_temperature)

        # apply deltas and equilibrate
        guess_temperature += delta_temperature
        guess_velocity += delta_velocity
        properties.set_equilibrium(
            working_gas=working_gas,
            density=guess_density,
            temperature=guess_temperature,
        )

    return guess_velocity


def _calc_density_ratio_speed(
    idx: int,
    current_density_ratio: float,
    initial_temperature: float,
    initial_pressure: float,
    species: SpeciesDefinition,
    mechanism: str,
    error_tol_temperature: float,
    error_tol_velocity: float,
):
    initial_state_gas = ct.Solution(mechanism)
    initial_state_gas.TPX = [
        initial_temperature,
        initial_pressure,
        species,
    ]
    working_gas = ct.Solution(mechanism)
    working_gas.TPX = [
        initial_temperature,
        initial_pressure,
        species,
    ]
    current_velocity = _calc_cj_guess_at_density_ratio_and_set_state(
        working_gas=working_gas,
        initial_state_gas=initial_state_gas,
        error_tol_temperature=error_tol_temperature,
        error_tol_velocity=error_tol_velocity,
        density_ratio=current_density_ratio,
    )

    return idx, current_velocity


def _calc_speed_parallel(
    initial_temperature: float,
    initial_pressure: float,
    species: SpeciesDefinition,
    mechanism: str,
    error_tol_temperature: float,
    error_tol_velocity: float,
    density_ratio_array: FloatArray,
) -> list[tuple[int, float]]:
    args = [
        [
            idx,
            ratio,
            initial_temperature,
            initial_pressure,
            species,
            mechanism,
            error_tol_temperature,
            error_tol_velocity,
        ]
        for idx, ratio in enumerate(density_ratio_array)
    ]
    with mp.Pool() as p:
        result = p.starmap(_calc_density_ratio_speed, args)

    return result


def _calc_speed_serial(
    initial_temperature: float,
    initial_pressure: float,
    species: SpeciesDefinition,
    mechanism: str,
    error_tol_temperature: float,
    error_tol_velocity: float,
    density_ratio_array: FloatArray,
) -> list[tuple[int, float]]:
    return list(
        map(
            _calc_density_ratio_speed,
            [item for item in range(len(density_ratio_array))],
            density_ratio_array,
            [initial_temperature for _ in density_ratio_array],
            [initial_pressure for _ in density_ratio_array],
            [species for _ in density_ratio_array],
            [mechanism for _ in density_ratio_array],
            [error_tol_temperature for _ in density_ratio_array],
            [error_tol_velocity for _ in density_ratio_array],
        )
    )


def shock(
    initial_pressure: float,
    initial_temperature: float,
    species: SpeciesDefinition,
    mechanism: str,
    parallelize=False,
) -> CjShock:
    """
    Calculates the Chapman-Jouguet detonation velocity of a gaseous mixture.

    This is a modified version of ``sdtoolbox.PostShock.CJspeed``. Specifically, it has been changed to allow for
    CJ state output, and to allow the use of multiprocessing to speed up the curve fit.

    :param initial_pressure: Initial pressure (Pa)
    :param initial_temperature: Initial temperature (K)
    :param species: Reactant species mole fractions
    :param mechanism: Cti file containing mechanism data, e.g. ``gri30.yaml``
    :param parallelize: Use multiprocessing for CJ state calculation, which is faster but requires the function to
        be run from ``__main__``, which may not behave well when used via Jupyter
    :return: Chapman-Jouguet detonation velocity (m/s)
    """
    num_steps = 20
    max_density_ratio = 2.0
    min_density_ratio = 1.5
    fit: CurveFit = CurveFit.empty()

    # Set error tolerances for CJ state calculation
    error_tol_temperature = 1e-4
    error_tol_velocity = 1e-4

    counter = 1
    r_squared = 0.0
    delta_r_squared = 0.0
    adjusted_density_ratio = 0.0

    while (counter <= 4) and ((r_squared < 0.99999) or (delta_r_squared < 1e-7)):
        density_ratios = np.linspace(min_density_ratio, max_density_ratio, num_steps + 1)

        if parallelize:
            speeds = _calc_speed_parallel(
                initial_temperature=initial_temperature,
                initial_pressure=initial_pressure,
                species=species,
                mechanism=mechanism,
                error_tol_temperature=error_tol_temperature,
                error_tol_velocity=error_tol_velocity,
                density_ratio_array=density_ratios,
            )

        else:
            speeds = _calc_speed_serial(
                initial_temperature=initial_temperature,
                initial_pressure=initial_pressure,
                species=species,
                mechanism=mechanism,
                error_tol_temperature=error_tol_temperature,
                error_tol_velocity=error_tol_velocity,
                density_ratio_array=density_ratios,
            )

        speeds.sort()
        cj_velocity_calculations = np.array(tuple(s for (_, s) in speeds))

        fit = curve_fit(density_ratios, cj_velocity_calculations)
        adjusted_density_ratio = -fit.b / (2.0 * fit.a)

        min_density_ratio = adjusted_density_ratio * (1 - 0.001)
        max_density_ratio = adjusted_density_ratio * (1 + 0.001)
        counter += 1

    cj_speed = fit.a * adjusted_density_ratio**2 + fit.b * adjusted_density_ratio + fit.c

    # calculate CJ shock state
    initial_state_gas = ct.Solution(mechanism)
    cj_state = ct.Solution(mechanism)
    initial_state_gas.TPX = [initial_temperature, initial_pressure, species]
    cj_state.TPX = [initial_temperature, initial_pressure, species]

    # All we want from this call is to mutate the working gas into CJ state; guess velocity is not needed.
    _ = _calc_cj_guess_at_density_ratio_and_set_state(
        working_gas=cj_state,
        initial_state_gas=initial_state_gas,
        error_tol_temperature=error_tol_temperature,
        error_tol_velocity=error_tol_velocity,
        density_ratio=adjusted_density_ratio,
    )

    return CjShock(speed=cj_speed, speed_r2=r_squared, state=cj_state)
