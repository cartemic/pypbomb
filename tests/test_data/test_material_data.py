from pypbomb.data import material
from pypbomb.data.material import MaterialProperties


def test_available_grades():
    assert "304L" in material.available_grades()


def test_all_material_groups():
    assert "2.1" in material.all_material_groups()


def test_properties():
    expected = MaterialProperties(
        group="2.3",
        elastic_modulus=197_500_000_000,
        density=7_970,
        poisson=0.27,
    )
    test = material.properties("316L")
    assert test == expected, "Faulty property read-in"
