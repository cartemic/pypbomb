import sqlite3
import typing
import warnings
from enum import Enum
from typing import Optional

import numpy as np

from . import _db, _units


class PipingType(Enum):
    seamless = "seamless"
    welded = "welded"


def max_stress_opt(piping_type: PipingType, grade: str, temp: float) -> Optional[float]:
    temp_f = _units.k_to_f(temp)
    if np.isnan(temp_f):
        warnings.warn(f"Invalid temperature: {temp} K", RuntimeWarning)
        return None

    with _db.connection() as conn:
        conn.row_factory = sqlite3.Row
        lower = conn.execute(
            """
            select
                temperature_f,
                max_stress_kpsi
            from b31_pipe_stress_limits
            where
                pipe_type=:pipe_type
                AND grade=:grade
                AND temperature_f <= :query_temp_f
            order by temperature_f desc
            limit 1;
            """,
            {
                "pipe_type": piping_type.value,
                "grade": grade,
                "query_temp_f": temp_f,
            },
        ).fetchone()
        upper = conn.execute(
            """
            select
                temperature_f,
                max_stress_kpsi
            from b31_pipe_stress_limits
            where
                pipe_type=:pipe_type
                AND grade=:grade
                AND temperature_f >= :query_temp_f
            order by temperature_f asc
            limit 1;
            """,
            {
                "pipe_type": piping_type.value,
                "grade": grade,
                "query_temp_f": temp_f,
            },
        ).fetchone()
        if None in (upper, lower):
            warnings.warn(
                f"No stress limits found for {piping_type.value} {grade} at {temp} K",
                RuntimeWarning,
            )
            return None

    temp_limits = [lower["temperature_f"], upper["temperature_f"]]
    stress_limits = [lower["max_stress_kpsi"], upper["max_stress_kpsi"]]

    return _units.kpsi_to_pa(typing.cast(float, np.interp(x=temp_f, xp=temp_limits, fp=stress_limits)))


def max_stress(piping_type: PipingType, grade: str, temp: float) -> float:
    maybe_stress = max_stress_opt(piping_type, grade, temp)
    if maybe_stress is None:
        raise ValueError(f"Unable to find max stress for piping type {piping_type.value}, grade {grade}, temp {temp}")
    return maybe_stress
