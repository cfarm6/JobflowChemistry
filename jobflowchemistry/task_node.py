from pydantic import dataclasses
from inspect import getfullargspec
from importlib import import_module
from enum import Enum
from monty.json import MontyDecoder, MontyEncoder
import json

class TaskNode():

    def as_dict(self) -> dict:
        """
        A JSON serializable dict representation of an object.
        """
        d: dict[str, Any] = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }

        try:
            parent_module = self.__class__.__module__.split(".", maxsplit=1)[0]
            module_version = import_module(parent_module).__version__
            d["@version"] = str(module_version)
        except (AttributeError, ImportError):
            d["@version"] = None

        spec = getfullargspec(self.__class__.__init__)

        def recursive_as_dict(obj):
            if isinstance(obj, (list, tuple)):
                return [recursive_as_dict(it) for it in obj]
            if isinstance(obj, dict):
                return {kk: recursive_as_dict(vv) for kk, vv in obj.items()}
            if hasattr(obj, "as_dict"):
                return obj.as_dict()
            if dataclasses is not None and dataclasses.is_dataclass(obj):
                d = dataclasses.asdict(obj)
                d.update(
                    {
                        "@module": obj.__class__.__module__,
                        "@class": obj.__class__.__name__,
                    }
                )
                return d
            return obj

        for c in spec.args + spec.kwonlyargs:
            if c != "self":
                try:
                    a = getattr(self, c)
                except AttributeError:
                    try:
                        a = getattr(self, "_" + c)
                    except AttributeError:
                        raise NotImplementedError(
                            "Unable to automatically determine as_dict "
                            "format from class. MSONAble requires all "
                            "args to be present as either self.argname or "
                            "self._argname, and kwargs to be present under "
                            "a self.kwargs variable to automatically "
                            "determine the dict format. Alternatively, "
                            "you can implement both as_dict and from_dict."
                        )
                d[c] = recursive_as_dict(a)
        if hasattr(self, "kwargs"):
            d.update(**self.kwargs)
        if spec.varargs is not None and getattr(self, spec.varargs, None) is not None:
            d.update({spec.varargs: getattr(self, spec.varargs)})
        if hasattr(self, "_kwargs"):
            d.update(**self._kwargs)
        if isinstance(self, Enum):
            d.update({"value": self.value})
        return d

    @classmethod
    def from_dict(cls, d):
        """

        Args:
            d: Dict representation.

        Returns:
            MSONable class.
        """
        decoded = {
            k: MontyDecoder().process_decoded(v)
            for k, v in d.items()
            if not k.startswith("@")
        }
        return cls(**decoded)

    def to_json(self) -> str:
        """
        Returns a json string representation of the MSONable object.
        """
        return json.dumps(self, cls=MontyEncoder)
