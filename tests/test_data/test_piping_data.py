from pypbomb.data import _units, piping


def test_available_schedules():
    good = {
        "8": {
            "5s",
            "5",
            "10s",
            "10",
            "20",
            "30",
            "40s",
            "40",
            "60",
            "80s",
            "80",
            "100",
            "120",
            "140",
            "160",
            "XXH",
        },
        "26": {"10", "20", "40s", "80s"},
    }
    for size, expected in good.items():
        test = piping.available_schedules(size)
        assert test == expected, f"Pipe schedule mismatch: expected {expected}, got {test}"


def test_dimensions():
    size = "9"
    good = {
        "40s": piping.PipeDimensions(
            inner_diameter=_units.in_to_m(8.941),
            outer_diameter=_units.in_to_m(9.625),
            mean_diameter=_units.in_to_m((8.941 + 9.625) / 2),
            wall_thickness=_units.in_to_m(0.342),
        ),
        "80s": piping.PipeDimensions(
            inner_diameter=_units.in_to_m(8.625),
            outer_diameter=_units.in_to_m(9.625),
            mean_diameter=_units.in_to_m((8.625 + 9.625) / 2),
            wall_thickness=_units.in_to_m(0.5),
        ),
    }
    for schedule, expected in good.items():
        test = piping.dimensions(size, schedule)
        assert test == expected, f"Pipe dimension mismatch: expected {expected}, got {test}"
