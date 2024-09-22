def temperature(t: float):
    if t < 0:
        raise ValueError(f"Invalid temperature: {t}. Must be > 0.")


def pressure(p: float):
    if p < 0:
        raise ValueError(f"Invalid pressure: {p}. Must be > 0.")
