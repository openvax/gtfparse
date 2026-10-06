"""Load the optional Polars conversion dependencies only when requested."""

import importlib


def require_polars():
    try:
        polars = importlib.import_module("polars")
    except ModuleNotFoundError as error:
        if error.name != "polars":
            raise
        raise ImportError(
            "Polars output requires the optional dependencies; "
            "install them with pip install 'gtfparse[polars]'"
        ) from error
    return polars
