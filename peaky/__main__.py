"""`python -m peaky ...` -> the same CLI as the `mascope-assign` script."""
import sys

from peaky.cli import main

if __name__ == "__main__":
    sys.exit(main())
