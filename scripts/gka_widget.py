"""Moved into the package -> peaky/gka_widget.py.

Preferred: `mascope-assign gka LEDGER.csv [-o out.html] [--ppm 2]`
       or  `python3 -m peaky.gka_widget LEDGER.csv ...`
This shim keeps the old `python3 scripts/gka_widget.py ...` invocation working.
"""
import os
import sys

# allow running straight from a checkout (no install needed)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from peaky.gka_widget import main  # noqa: E402

if __name__ == "__main__":
    main()
