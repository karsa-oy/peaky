"""Back-compat shim: ``mascope_assign`` is now ``peaky``.

The package was renamed to ``peaky`` (its distribution + CLI name). This module
keeps the old import path working as a true alias of the real package — including
submodules — so existing ``import mascope_assign`` / ``from mascope_assign import
x`` / ``import mascope_assign.cli`` code is unchanged AND resolves to the SAME
module objects as ``peaky`` (no duplicate, separately-stateful module copies).

How: a meta-path finder redirects any ``mascope_assign[.sub]`` import to the
corresponding ``peaky[.sub]`` module, then ``sys.modules`` is pointed at ``peaky``
so bare ``import mascope_assign`` and attribute access also hit the real package.
Keep this file a pure alias — ``peaky`` is the single source of truth.
"""
import importlib
import importlib.abc
import importlib.util
import sys

_OLD = __name__          # "mascope_assign"
_NEW = "peaky"


class _Redirector(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Resolve ``mascope_assign`` and ``mascope_assign.<sub>`` to the live
    ``peaky`` / ``peaky.<sub>`` module object."""

    def find_spec(self, fullname, path=None, target=None):
        if fullname == _OLD or fullname.startswith(_OLD + "."):
            return importlib.util.spec_from_loader(fullname, self)
        return None

    def create_module(self, spec):
        new_name = _NEW + spec.name[len(_OLD):]   # mascope_assign.cli -> peaky.cli
        module = importlib.import_module(new_name)
        sys.modules[spec.name] = module           # alias points at the real object
        return module

    def exec_module(self, module):                # already executed as peaky.*
        pass


if not any(isinstance(f, _Redirector) for f in sys.meta_path):
    sys.meta_path.insert(0, _Redirector())

import peaky  # noqa: E402

sys.modules[__name__] = peaky
