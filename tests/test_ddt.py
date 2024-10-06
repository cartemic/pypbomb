import cantera as ct
import numpy as np
import pytest

from pypbomb import ddt


def test_calculate_spiral_diameter():
    expected = 0.0184100557685425  # hand calc, meters
    result = ddt.ideal_blockage_diameter(
        pipe_id=0.146304,  # 5.76 in -> m
        br=0.44,
    )

    assert np.isclose(result, expected)


def test_calculate_ddt_run_up_bad_blockage_ratio():
    bad_blockages = [-4.0, 0, 1]
    for br in bad_blockages:
        with pytest.raises(ValueError, match="Blockage ratio outside of correlation range"):
            ddt.runup_distance(
                pipe_id=1,
                br=br,
                initial_temperature=300,
                initial_pressure=101325,
                species={"H2": 1, "O2": 0.5},
                mechanism="gri30.cti",
                phase_specification="gri30",
            )


def test_calculate_ddt_run_up():
    mechanism = "gri30.yaml"
    gas = ct.Solution(mechanism)
    initial_temperature = 300  # K
    initial_pressure = 101_325  # Pa
    gas.set_equivalence_ratio(1, "CH4", {"O2": 1, "N2": 3.76})
    species = gas.mole_fraction_dict()
    # expected numbers approximated from fig 45
    conditions = {
        "eq4.1, upper bound": (0.1, 48.51385390428211),
        "interpolated": (0.2, 29.24433249370277),
        "eq4.4, lower bound": (0.3, 18.136020151133494),
        "eq4.4, upper bound": (0.75, 4.76070528967254),
    }

    for test_condition, (br, expected) in conditions.items():
        result = ddt.runup_distance(
            pipe_id=1,  # figure specifies 1m diameter, so Xs/D is equivalent in magnitude to Xs
            br=br,
            initial_temperature=initial_temperature,
            initial_pressure=initial_pressure,
            species=species,
            mechanism=mechanism,
            phase_specification="gri30",
        )
        # Tolerance on this test is a bit loose because the x-axis of fig 45 is logarithmic.
        # Additionally, we want a little wiggle room to account for differences in simulations (mechanism, etc.).
        assert np.isclose(expected, result, rtol=0.5)
