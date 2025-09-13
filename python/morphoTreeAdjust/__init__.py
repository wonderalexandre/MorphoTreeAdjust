try:
    from .version import __version__  # opcional
except Exception:  # pragma: no cover
    __version__ = "0.0.0"

from .morphoTreeAdjust import *  # tipos e classes C++ (pybind)

# Helpers Python movidos para submódulo: `from morphoTreeAdjust import utils as mta_utils`
# Obs.: `viz` permanece como alias compatível.

__all__ = [
    '__version__',
]
