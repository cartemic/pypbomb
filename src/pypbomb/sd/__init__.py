"""
Pypbomb requires the use of some functionality from the shock and detonation toolbox from the Explosion Dynamics
Laboratory at Cal Tech [1]. SDToolbox is, as yet, not installable via pip or conda, and not all of its functions are
needed for this package. Therefore, the functions needed for pypbomb's purposes are located here, in ``pypbomb.sd``.
Some of these functions, such as ``Detonation.cj_speed``, have been modified to speed up the calculation process.
If you are planning to do any sort of actual detonation work outside of tube design, it is recommended to use the full
SDToolbox.

[1] SDToolbox Numerical Tools for Shock and Detonation Wave Modeling,
    Explosion Dynamics Laboratory, Contributors: S. Browne, J. Ziegler,
    N. Bitter, B. Schmidt, J. Lawson and J. E. Shepherd
    California Institute of Technology, Pasadena, CA USA
    http://shepherd.caltech.edu/EDL/PublicResources/sdt/
"""

from . import cj, error, properties, reflect
