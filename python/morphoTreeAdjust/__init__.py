"""Python entrypoint for the MorphoTreeAdjust core API."""

try:
    from .version import __version__  # optional
except Exception:  # pragma: no cover
    __version__ = "0.0.0"

from .morphoTreeAdjust import *  # C++ types and classes (pybind)

# Python helpers moved to a submodule: `from morphoTreeAdjust import utils as mta_utils`

__all__ = [
    '__version__',
    'AdjacencyRelation',
    'DynamicComponentTree',
    'DualMinMaxTreeIncrementalFilter',
    'ComponentTreeCasf',
]
