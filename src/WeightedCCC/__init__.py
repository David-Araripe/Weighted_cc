# -*- coding: utf-8 -*-
try:
    from .graph import GraphClosure  # noqa: F401
except ImportError:
    pass

from ._version_helper import get_version

__version__ = get_version()
