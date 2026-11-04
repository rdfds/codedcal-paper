"""Allow ``python -m codedcal`` as an alias for the ``codedcal`` script."""
from .cli import main

if __name__ == '__main__':
    raise SystemExit(main())
