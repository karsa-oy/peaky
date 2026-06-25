"""The assignment-pass director, split into cohesive submodules.

Public surface is unchanged: every name the pipeline used as `passes.X`
is re-exported here (config -> core -> postprocess -> directors)."""

from .config import *  # noqa: F401,F403
from .core import *  # noqa: F401,F403
from .postprocess import *  # noqa: F401,F403
from .directors import *  # noqa: F401,F403

__version__ = "0.8.0"  # offset-tolerant gates + prior_offset + calibration-aware arbitration + carbon-clamp Si-skip
