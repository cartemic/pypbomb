import numpy as np

from pypbomb.data import _units


def test_in_to_m():
    assert np.isnan(_units.in_to_m(-8)), "nonphysical result is nan (scalar)"
    assert np.isclose(_units.in_to_m(1), 0.0254), "correct conversion (scalar)"
    assert np.allclose(
        _units.in_to_m(np.array([1, -8])),
        [0.0254, np.nan],
        equal_nan=True,
    ), "nonphysical result is nan (array)"
    assert np.allclose(_units.in_to_m(np.array([1, 1])), [0.0254, 0.0254]), "correct conversion (array)"


def test_tpi_to_tpm():
    assert np.isnan(_units.tpi_to_tpm(-8)), "nonphysical result is nan (scalar)"
    assert np.isclose(_units.tpi_to_tpm(1), 39.3700787401), "correct conversion (scalar)"
    assert np.allclose(
        _units.tpi_to_tpm(np.array([1, -8])),
        [39.3700787401, np.nan],
        equal_nan=True,
    ), "nonphysical result is nan (array)"
    assert np.allclose(
        _units.tpi_to_tpm(np.array([1, 1])),
        [39.3700787401, 39.3700787401],
    ), "correct conversion (array)"


def test_gpa_to_pa():
    assert np.isnan(_units.gpa_to_pa(-8)), "nonphysical result is nan (scalar)"
    assert np.isclose(_units.gpa_to_pa(17), 17_000_000_000), "correct conversion (scalar)"
    assert np.allclose(
        _units.gpa_to_pa(np.array([17, -8])),
        [17_000_000_000, np.nan],
        equal_nan=True,
    ), "nonphysical result is nan (array)"
    assert np.allclose(
        _units.gpa_to_pa(np.array([17, 17])),
        [17_000_000_000, 17_000_000_000],
    ), "correct conversion (array)"


def test_kpsi_to_pa():
    assert np.isnan(_units.kpsi_to_pa(-8)), "nonphysical result is nan (scalar)"
    assert np.isclose(_units.kpsi_to_pa(12), 82_737_087.518), "correct conversion (scalar)"
    assert np.allclose(
        _units.kpsi_to_pa(np.array([12, -8])),
        [82_737_087.518, np.nan],
        equal_nan=True,
    ), "nonphysical result is nan (array)"
    assert np.allclose(
        _units.kpsi_to_pa(np.array([12, 12])),
        [82_737_087.518, 82_737_087.518],
    ), "correct conversion (array)"


def test_psi_to_pa():
    assert np.isnan(_units.psi_to_pa(-8)), "nonphysical result is nan (scalar)"
    assert np.isclose(_units.psi_to_pa(12), 82_737.087518), "correct conversion (scalar)"
    assert np.allclose(
        _units.psi_to_pa(np.array([12, -8])),
        [82_737.087518, np.nan],
        equal_nan=True,
    ), "nonphysical result is nan (array)"
    assert np.allclose(
        _units.psi_to_pa(np.array([12, 12])),
        [82_737.087518, 82_737.087518],
    ), "correct conversion (array)"


def test_gcm3_to_kgm3():
    assert np.isnan(_units.gcm3_to_kgm3(-8)), "nonphysical result is nan (scalar)"
    assert np.isclose(_units.gcm3_to_kgm3(8), 8_000), "correct conversion (scalar)"
    assert np.allclose(
        _units.gcm3_to_kgm3(np.array([8, -8])), [8_000, np.nan], equal_nan=True
    ), "nonphysical result is nan (array)"
    assert np.allclose(_units.gcm3_to_kgm3(np.array([8, 8])), [8_000, 8_000]), "correct conversion (array)"


def test_k_to_f():
    assert np.isnan(_units.k_to_f(-8)), "nonphysical result is nan (scalar)"
    assert np.isclose(_units.k_to_f(19), -425.47), "correct conversion (scalar)"
    assert np.allclose(
        _units.k_to_f(np.array([19, -8])),
        [-425.47, np.nan],
        equal_nan=True,
    ), "nonphysical result is nan (array)"
    assert np.allclose(
        _units.k_to_f(np.array([19, 19])),
        [-425.47, -425.47],
    ), "correct conversion (array)"


def test_f_to_k():
    assert np.isnan(_units.f_to_k(-500)), "nonphysical result is nan (scalar)"
    assert np.isclose(_units.f_to_k(70), 294.2611111111112), "correct conversion (scalar)"
    assert np.allclose(
        _units.f_to_k(np.array([70, -500])),
        [294.2611111111112, np.nan],
        equal_nan=True,
    ), "nonphysical result is nan (array)"
    assert np.allclose(
        _units.f_to_k(np.array([70, 70])),
        [294.2611111111112, 294.2611111111112],
    ), "correct conversion (array)"
