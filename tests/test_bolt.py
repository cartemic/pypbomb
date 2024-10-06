import numpy as np

from pypbomb import bolt
from pypbomb.bolt import BoltFailureMode, BoltStressAreas, StressAreas
from pypbomb.data.bolt import (
    DiameterRange,
    ExternalDiameterAnsi,
    ExternalThreadAnsi,
    InternalDiameterAnsi,
    InternalThreadAnsi,
)


def test_safety_factors():
    result = bolt.safety_factors(
        internal_pressure=11250175.475279031,  # 1631.7 psi -> Pa,
        external_pressure=101_325,
        window_area=0.009274175,  # 14.375 in^2 (5.75x2.5) -> m^2
        num_bolts=20,
        thread_size="1/4-28",
        partial_thread_class="2",
        bolt_max_tensile=827_370_875.1802038,  # 120 kpsi -> Pa
        plate_max_tensile=206_842_718.79505095,  # 30 kpsi -> Pa
        engagement_length=0.0127,  # 0.5 in -> m
    )

    assert np.isclose(result.bolt, 3.606968028), "bolt"
    assert np.isclose(result.plate, 7.969517321), "plate"


def test_bolt_area_shear():
    tpm = 1_102.36220472  # 28 tpi -> tpm
    engagement_length = 0.0127  # 0.5 in -> m
    k_n_max = 0.005588  # 0.22 in -> m
    e_s_min = 0.0056515  # 0.2225 in -> m
    expected = 0.00012048611034360326  # m^2
    result = bolt._bolt_area_shear(tpm=tpm, engagement_length=engagement_length, k_n_max=k_n_max, e_s_min=e_s_min)

    assert np.isclose(result, expected)


def test_bolt_area_tensile_lt_100kpsi():
    tpm = 1_102.36220472  # 28 tpi -> tpm
    d_basic = 0.00635  # 0.25 in -> m
    expected = 2.3467096683689698e-5  # m^2
    result = bolt._bolt_area_tensile_lt_100kpsi(tpm=tpm, d_basic=d_basic)

    assert np.isclose(result, expected)


def test_bolt_area_tensile_ge_100kpsi():
    tpm = 1_102.36220472  # 28 tpi -> tpm
    e_s_min = 0.0056515  # 0.2225 in -> m
    expected = 2.2538051016840217e-5  # m^2
    result = bolt._bolt_area_tensile_ge_100kpsi(tpm=tpm, e_s_min=e_s_min)

    assert np.isclose(result, expected)


def test_plate_area_shear():
    tpm = 1_102.36220472  # 28 tpi -> tpm
    engagement_length = 0.0127  # 0.5 in -> m
    d_s_min = 0.0063246  # 0.249 in -> m
    e_n_max = 0.00586994  # 0.2311 in -> m
    expected = 0.00019918933192013787  # m^2
    result = bolt._plate_area_shear(tpm=tpm, engagement_length=engagement_length, d_s_min=d_s_min, e_n_max=e_n_max)

    assert np.isclose(result, expected)


def test_stress_areas():
    result = bolt.stress_areas(
        thread_size="1/4-28",
        partial_thread_class="2",
        bolt_max_tensile=827_370_875.1802038,  # 120 kpsi -> Pa
        engagement_length=0.0127,  # 0.5 in -> m
    )

    assert np.isclose(result.plate, 0.00019918933192013787)
    assert np.isclose(result.bolt.shear, 0.00012048611034360326)
    assert np.isclose(result.bolt.tensile, 2.2538051016840217e-5)
    assert result.bolt.failure_mode is BoltFailureMode.tension


def test_minimum_engagement():
    bolt_max_tensile = 827_370_875.1802038  # 120 kpsi -> Pa
    plate_max_tensile = 206_842_718.79505095  # 30 kpsi -> Pa
    bolt_dimensions = ExternalThreadAnsi(
        thread_size="1/4-28",
        thread_class="2A",
        threads_per_inch=28,
        diameter=ExternalDiameterAnsi(
            basic=0.00635,
            major=DiameterRange(min=0.0063246, max=0.006159499999999999),
            pitch=DiameterRange(min=0.0056515, max=0.0057353199999999995),
            minor=0.00521208,
        ),
    )
    plate_dimensions = InternalThreadAnsi(
        thread_size="1/4-28",
        thread_class="2B",
        threads_per_inch=28,
        diameter=InternalDiameterAnsi(
            basic=0.00635,
            minor=DiameterRange(min=0.005359399999999999, max=0.005588),
            pitch=DiameterRange(min=0.0057607199999999996, max=0.00586994),
            major=0.00635,
        ),
    )
    stress_areas = StressAreas(
        bolt=BoltStressAreas(
            shear=0.00012048611034360326,  # m^2
            tensile=2.2538051016840217e-5,  # m^2
            failure_mode=BoltFailureMode.tension,
        ),
        plate=0.00019918933192013787,  # m^2
    )
    expected = 0.011495926821166583  # m
    result = bolt.minimum_engagement(
        bolt_dimensions=bolt_dimensions,
        plate_dimensions=plate_dimensions,
        areas=stress_areas,
        bolt_max_tensile=bolt_max_tensile,
        plate_max_tensile=plate_max_tensile,
    )

    assert np.isclose(result, expected)


def test_j_factor():
    bolt_max_tensile = 827_370_875.1802038  # 120 kpsi -> Pa
    plate_max_tensile = 206_842_718.79505095  # 30 kpsi -> Pa
    stress_areas = StressAreas(
        bolt=BoltStressAreas(
            shear=0.00012048611034360326,  # m^2
            tensile=2.2538051016840217e-5,  # m^2
            failure_mode=BoltFailureMode.tension,
        ),
        plate=0.00019918933192013787,  # m^2
    )
    expected = 2.419529383067773  # unitless
    result = bolt._j_factor(
        areas=stress_areas,
        bolt_max_tensile=bolt_max_tensile,
        plate_max_tensile=plate_max_tensile,
    )

    assert np.isclose(result, expected)
