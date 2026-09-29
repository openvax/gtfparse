import os


def data_path(name):
    """
    Return the absolute path to a fixture relative to tests/data.
    """
    return os.path.join(os.path.dirname(__file__), "data", name)
