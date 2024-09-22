from typing import Any

import numpy as np

FloatArray = np.ndarray[Any, np.dtype[np.floating]]
MoleFractions = str | dict[str, float]
