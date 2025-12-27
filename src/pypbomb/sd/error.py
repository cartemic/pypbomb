from dataclasses import dataclass

from cantera import Solution


@dataclass(frozen=True)
class StateError:
    pressure: float
    enthalpy: float


def equilibrium(
    working_gas: Solution,
    initial_state_gas: Solution,
    speed: float,
) -> StateError:
    """
    Calculate error in pressure and enthalpy using the momentum and energy conservation equations for equilibrium
    calculations. Nothing is mutated.

    :param working_gas: Working gas mixture used for calculations
    :param initial_state_gas: Working gas mixture in its initial, undetonated state
    :param speed: A guess for the initial velocity (m/s)
    :return:
    """
    initial_pressure = initial_state_gas.P
    initial_enthalpy = initial_state_gas.enthalpy_mass
    initial_density = initial_state_gas.density
    initial_velocity = speed

    working_pressure = working_gas.P
    working_enthalpy = working_gas.enthalpy_mass
    working_density = working_gas.density

    working_velocity = initial_velocity * initial_density / working_density

    sqr_vel_initial = initial_velocity**2
    sqr_vel_working = working_velocity**2

    enthalpy_error = (working_enthalpy + 0.5 * sqr_vel_working) - (initial_enthalpy + 0.5 * sqr_vel_initial)
    pressure_error = (working_pressure + working_density * sqr_vel_working) - (
        initial_pressure + initial_density * sqr_vel_initial
    )

    return StateError(pressure=pressure_error, enthalpy=enthalpy_error)


def reflected_shock_frozen(
    working_gas: Solution,
    post_shock_gas: Solution,
    speed: float,
) -> StateError:
    """
    Calculate error in pressure and enthalpy using the momentum and energy conservation equations for reflected shock
    calculations. Nothing is mutated.

    :param working_gas: Working gas mixture used for calculations
    :param post_shock_gas: Working gas mixture in its post-incident-shock state
    :param speed: Current post-incident-shock lab frame particle speed (m/s)
    """
    post_shock_pressure = post_shock_gas.P
    post_shock_enthalpy = post_shock_gas.enthalpy_mass
    post_shock_density = post_shock_gas.density

    working_pressure = working_gas.P
    working_enthalpy = working_gas.enthalpy_mass
    working_density = working_gas.density

    enthalpy_error = (
        working_enthalpy
        - post_shock_enthalpy
        - 0.5 * (speed**2) * ((working_density / post_shock_density) + 1) / (working_density / post_shock_density - 1)
    )
    pressure_error = (
        working_pressure
        - post_shock_pressure
        - working_density * (speed**2) / (working_density / post_shock_density - 1)
    )

    return StateError(pressure=pressure_error, enthalpy=enthalpy_error)
