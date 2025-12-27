import sqlite3
from dataclasses import dataclass
from typing import Optional

from . import _db, _units


@dataclass
class MaterialProperties:
    group: str
    """
    ASME B16.34 material group
    """
    elastic_modulus: float
    """
    Elastic modulus (Pa)
    """
    density: float
    """
    Mass density (kg/m^3)
    """
    poisson: float
    """
    Poisson ratio (unitless)
    """


def available_grades() -> set[str]:
    """
    :return: Available material grades
    """
    with _db.connection() as conn:
        return set(result[0] for result in conn.execute("select distinct grade from material;").fetchall())


def all_material_groups() -> set[str]:
    """
    :return: Available ASME B16.34 material groups
    """
    with _db.connection() as conn:
        return set(result[0] for result in conn.execute("select distinct material_group from material;").fetchall())


def properties_opt(grade: str) -> Optional[MaterialProperties]:
    """
    :param grade: Material grade, e.g. "316L"
    :return: Material properties
    """
    with _db.connection() as conn:
        conn.row_factory = sqlite3.Row
        result = conn.execute(
            """
            select
                material_group,
                elastic_modulus_gpa,
                density_g_cm3,
                poisson
            from material
            where grade=:grade;
            """,
            {"grade": grade},
        ).fetchone()
    if result is None:
        return result
    return MaterialProperties(
        group=result["material_group"],
        elastic_modulus=_units.gpa_to_pa(result["elastic_modulus_gpa"]),
        density=_units.gcm3_to_kgm3(result["density_g_cm3"]),
        poisson=result["poisson"],
    )


def properties(grade: str) -> MaterialProperties:
    """
    :param grade: Material grade, e.g. "316L"
    :return: Material properties
    """
    maybe_properties = properties_opt(grade)
    if maybe_properties is None:
        raise ValueError(f"No material properties found for grade {grade}")
    return maybe_properties
