"""Deprecated entry point — superseded by the installed console command:

    peaky assign --sample-id <ID> --reagent <Br|Ur|auto> ...
    (or: python3 -m peaky assign ...)

Kept as a thin forwarder so the old `python3 scripts/run_assignment.py ...` call
still works. The flat flags are forwarded verbatim to the `assign` subcommand,
which now also accepts --reagent / --adducts / --env. Heavy work runs on the host
Python (this package + mascope-sdk); run via the shell.
"""
import os
import sys

# allow running straight from a checkout (no install needed)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from peaky.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main(["assign", *sys.argv[1:]]))
