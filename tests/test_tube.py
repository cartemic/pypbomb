import numpy as np

from pypbomb import thermochem, tube
from pypbomb.data.material import MaterialProperties
from pypbomb.data.piping import PipeDimensions
from pypbomb.tube import TubeKind


def test_dynamic_load_factor():
    # hand calcs use:
    #     material = "316L"
    #     schedule = "80"
    #     nominal_size = "6"
    tube_id = 0.1463294  # 5.761 in -> m
    tube_od = 0.168275  # 6.625 in -> m
    pipe_dimensions = PipeDimensions(
        inner_diameter=tube_id,
        outer_diameter=tube_od,
        mean_diameter=(tube_od + tube_id) / 2,
        wall_thickness=(tube_od - tube_id) / 2,
    )
    material_properties = MaterialProperties(
        group="2.3",
        elastic_modulus=200_000_000_000,  # 200 GPa -> Pa
        density=7_900,  # 7.9 g/cm**3 -> kg/m**3
        poisson=0.28,
    )

    # from hand calcs, critical velocity is 1457.44 m/s, giving upper and
    # lower bounds of 1603.188 and 1311.700 m/s
    test_values = {
        # CJ speed: DLF
        1200: 1,
        1311: 1,
        1312: 4,
        1400: 4,
        1603: 4,
        1604: 2,
        2000: 2,
    }
    for cj_velocity, expected in test_values.items():
        test = tube.dynamic_load_factor(
            pipe_dimensions=pipe_dimensions,
            material_properties=material_properties,
            cj_velocity=cj_velocity,
            critical_velocity_band_frac=0.1,
        )
        assert np.isclose(test, expected), f"Incorrect DLF calculation for CJ speed {cj_velocity}"


def test_calculate_max_pressure():
    pipe_dimensions = PipeDimensions(
        wall_thickness=0.0109728,  # 0.432 in -> m
        outer_diameter=0.168275,  # 6.625 in -> m
        inner_diameter=0.1463294,  # 5.761 in -> m
        mean_diameter=0.1573022,  # 6.193 in -> m
    )
    test = tube.max_allowable_pressure(
        pipe_dimensions=pipe_dimensions,
        max_stress=129621436.86400001,  # 18.8 kpsi -> Pa
        safety_factor=4,
    )
    assert np.isclose(test, 4520947.90289424)


def test_calculate_initial_pressure_closed_tube():
    # NPS 6 schedule 80
    initial_temperature = 300  # K
    tube_id = 0.1463294  # 5.761 in -> m
    tube_od = 0.168275  # 6.625 in -> m
    reactant_mixture = dict(H2=1 / 3, O2=2 / 3)
    mechanism = "gri30.yaml"
    rtol = 1e-4

    # the initial pressure should cause the reflected detonation pressure
    # to be equal to the tube's max pressure, accounting for dynamic load
    # factor
    max_pressure = 15101889.83212825  # Pa
    test = tube.max_initial_pressure(
        max_pressure=max_pressure,
        initial_temperature=300,
        mechanism=mechanism,
        species=reactant_mixture,
        pipe_dimensions=PipeDimensions(
            inner_diameter=tube_id,
            outer_diameter=tube_od,
            wall_thickness=(tube_od - tube_id) / 2,
            mean_diameter=(tube_od + tube_id) / 2,
        ),
        # 316L
        material_properties=MaterialProperties(
            group="2.2",
            elastic_modulus=200_000_000_000,  # 200 GPa -> Pa
            density=7_900,  # g/cm^3 -> kg/m^3
            poisson=0.28,
        ),
        parallelize_cj_calc=True,
        tube_kind=TubeKind.closed,
        rtol=rtol,
    )
    # Use calculated initial pressure for a reflected shock calculation and make sure we get the DLF-scaled maximum
    # pressure as a result
    reflected_pressure = thermochem.reflected_cj_shock(
        initial_temperature=initial_temperature,
        initial_pressure=test.pressure,
        species=reactant_mixture,
        mechanism=mechanism,
    ).reflected.state.P
    error = abs(max_pressure - reflected_pressure * test.dlf) / max_pressure
    assert error <= rtol


def test_calculate_initial_pressure_open_tube():
    # NPS 6 schedule 80
    initial_temperature = 300  # K
    tube_id = 0.1463294  # 5.761 in -> m
    tube_od = 0.168275  # 6.625 in -> m
    reactant_mixture = dict(H2=1 / 3, O2=2 / 3)
    mechanism = "gri30.yaml"
    rtol = 1e-4

    # the initial pressure should cause the reflected detonation pressure
    # to be equal to the tube's max pressure, accounting for dynamic load
    # factor
    max_pressure = 15101889.83212825  # Pa
    test = tube.max_initial_pressure(
        max_pressure=max_pressure,
        initial_temperature=300,
        mechanism=mechanism,
        species=reactant_mixture,
        pipe_dimensions=PipeDimensions(
            inner_diameter=tube_id,
            outer_diameter=tube_od,
            wall_thickness=(tube_od - tube_id) / 2,
            mean_diameter=(tube_od + tube_id) / 2,
        ),
        # 316L
        material_properties=MaterialProperties(
            group="2.2",
            elastic_modulus=200_000_000_000,  # 200 GPa -> Pa
            density=7_900,  # g/cm^3 -> kg/m^3
            poisson=0.28,
        ),
        parallelize_cj_calc=True,
        tube_kind=TubeKind.open,
        rtol=rtol,
    )
    # Use calculated initial pressure for a CJ shock calculation and make sure we get the DLF-scaled maximum
    # pressure as a result
    reflected_pressure = thermochem.reflected_cj_shock(
        initial_temperature=initial_temperature,
        initial_pressure=test.pressure,
        species=reactant_mixture,
        mechanism=mechanism,
    ).cj.state.P
    error = abs(max_pressure - reflected_pressure * test.dlf) / max_pressure
    assert error <= rtol
