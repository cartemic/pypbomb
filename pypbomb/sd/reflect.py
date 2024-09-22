import warnings
from dataclasses import dataclass

import cantera as ct
import numpy as np
from cantera import Solution

from pypbomb.sd import error, properties


@dataclass(frozen=True)
class ReflectedShock:
    speed: float
    state: ct.Solution


def shock(
    mechanism: str,
    initial_state_gas: Solution,
    post_shock_gas: Solution,
    incident_shock_speed: float,
) -> ReflectedShock:
    """
    Calculate equilibrium post-reflected-shock state assuming u1 = 0.

    Original function: ``sdtoolbox.reflections.reflected_eq``

    :param mechanism: Cti file containing mechanism data (e.g. ``gri30.yaml``)
    :param initial_state_gas: Working gas mixture in its initial, undetonated state.
    :param post_shock_gas: Working gas mixture in its post-incident-shock state.
    :param incident_shock_speed: Incident shock speed in m/s.
    :return:
    """
    initial_pressure = initial_state_gas.P
    initial_volume = 1 / initial_state_gas.density

    reflected_pressure = post_shock_gas.P
    reflected_density = post_shock_gas.density
    reflected_volume = 1 / post_shock_gas.density
    reflected_temperature = post_shock_gas.T
    reflected_velocity = np.sqrt((reflected_pressure - initial_pressure) * (initial_volume - reflected_volume))

    working_volume = 0.2 / reflected_density
    working_pressure = reflected_pressure + reflected_density * (incident_shock_speed**2) * (
        1 - working_volume / reflected_volume
    )
    working_temperature = (
        reflected_temperature * working_pressure * working_volume / (reflected_pressure * reflected_volume)
    )

    working_gas = ct.Solution(mechanism)
    working_gas.TPX = [working_temperature, working_pressure, post_shock_gas.X]
    _set_eq_state(reflected_velocity, post_shock_gas, working_gas)
    working_pressure = working_gas.P
    reflected_shock_speed = (
        working_pressure - reflected_pressure
    ) / reflected_velocity / reflected_density - reflected_velocity

    return ReflectedShock(state=working_gas, speed=reflected_shock_speed)


def _set_eq_state(
    particle_speed: float,
    post_shock_gas: ct.Solution,
    working_gas: ct.Solution,
    error_tol_temperature: float = 1e-4,
    error_tol_specific_volume: float = 1e-4,
    max_iterations: int = 500,
):
    """
    Set the ``working_gas`` object to the equilibrium post-reflected-shock state for a shock traveling at
    ``particle_speed``

    Original function: ``sdtoolbox.reflections.PostReflectedShock_eq``

    :param particle_speed: Post-incident-shock lab frame particle speed (m/s)
    :param post_shock_gas: Working gas mixture in its post-incident-shock state
    :param working_gas: Working gas mixture used for calculations (MUTATED)
    :param error_tol_temperature: Temperature error tolerance for iteration (K)
    :param error_tol_specific_volume: Specific volume error tolerance for iteration (m^3/kg)
    :param max_iterations: Maximum number of loop iterations
    :return:
    """
    post_shock_volume = 1 / post_shock_gas.density

    guess_temperature = working_gas.T
    guess_density = working_gas.density
    guess_volume = 1 / guess_density
    properties.set_equilibrium(working_gas=working_gas, density=guess_density, temperature=guess_temperature)

    delta_temperature = 1000
    delta_volume = 1000

    loop_counter = 0
    while (abs(delta_temperature) > error_tol_temperature * guess_temperature) or (
        abs(delta_volume) > error_tol_specific_volume * guess_volume
    ):
        loop_counter += 1
        if loop_counter == max_iterations:  # pragma: no cover
            warnings.warn(
                "Calculation did not converge for U = {0:.2f}after {1} iterations".format(particle_speed, loop_counter)
            )
            return working_gas

        # calculate enthalpy and pressure error for current guess
        err_unperturbed = error.reflected_shock_frozen(
            working_gas=working_gas,
            post_shock_gas=post_shock_gas,
            speed=particle_speed,
        )

        # equilibrate working gas with perturbed temperature
        delta_temperature = guess_temperature * 0.02
        properties.set_equilibrium(
            working_gas=working_gas,
            density=guess_density,
            temperature=guess_temperature + delta_temperature,
        )

        # calculate enthalpy and pressure error for perturbed temperature
        err_perturbed = error.reflected_shock_frozen(
            working_gas=working_gas,
            post_shock_gas=post_shock_gas,
            speed=particle_speed,
        )

        # calculate temperature derivatives
        deriv_enthalpy_temperature = (err_perturbed.enthalpy - err_unperturbed.enthalpy) / delta_temperature
        deriv_pressure_temperature = (err_perturbed.pressure - err_unperturbed.pressure) / delta_temperature

        # equilibrate working gas with perturbed volume
        delta_volume = 0.02 * guess_volume
        properties.set_equilibrium(
            working_gas=working_gas,
            density=1 / (guess_volume + delta_volume),
            temperature=guess_temperature,
        )

        # calculate enthalpy and pressure error for perturbed specific vol
        err_perturbed = error.reflected_shock_frozen(
            working_gas=working_gas,
            post_shock_gas=post_shock_gas,
            speed=particle_speed,
        )

        # calculate specific volume derivatives
        deriv_enthalpy_volume = (err_perturbed.enthalpy - err_unperturbed.enthalpy) / delta_volume
        deriv_pressure_volume = (err_perturbed.pressure - err_unperturbed.pressure) / delta_volume

        # solve matrix for temperature and volume deltas
        j = deriv_enthalpy_temperature * deriv_pressure_volume - deriv_pressure_temperature * deriv_enthalpy_volume
        bb = [
            deriv_pressure_volume,
            -deriv_enthalpy_volume,
            -deriv_pressure_temperature,
            deriv_enthalpy_temperature,
        ]
        aa = [-err_unperturbed.enthalpy, -err_unperturbed.pressure]

        delta_temperature = (bb[0] * aa[0] + bb[1] * aa[1]) / j
        delta_volume = (bb[2] * aa[0] + bb[3] * aa[1]) / j

        # check and limit temperature delta
        delta_temp_max = 0.2 * guess_temperature
        if abs(delta_temperature) > delta_temp_max:
            delta_temperature = delta_temp_max * delta_temperature / abs(delta_temperature)

        # check and limit specific volume delta
        perturbed_volume = guess_volume + delta_volume
        if perturbed_volume > post_shock_volume:
            delta_volume_max = 0.5 * (post_shock_volume - guess_volume)
        else:
            delta_volume_max = 0.2 * guess_volume

        if abs(delta_volume) > delta_volume_max:
            delta_volume = delta_volume_max * delta_volume / abs(delta_volume)

        # apply calculated and limited deltas to temperature and spec. vol
        guess_temperature += +delta_temperature
        guess_volume += delta_volume
        guess_density = 1 / guess_volume

        properties.set_equilibrium(
            working_gas=working_gas,
            density=guess_density,
            temperature=guess_temperature,
        )
