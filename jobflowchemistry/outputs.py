from pydantic.dataclasses import dataclass

from .Structure import Structure

class Settings(dict):
    """
    Dictionary-like container for calculation settings.
    """
    pass

class Properties(dict):
    """
    Dictionary-like container for computed properties.
    """
    pass
