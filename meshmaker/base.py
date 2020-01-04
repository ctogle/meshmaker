class Base:
    """Base class readily accepts defaults"""

    def __init__(self, **kws):
        for k, v in kws.items():
            if not hasattr(self, k):
                setattr(self, k, v)


class Laziness(Base):
    """Generally lazy lookup"""

    def __init__(self, method, **kws):
        super().__init__(**kws)
        self._lookup = {}
        self._method = method

    def __getitem__(self, key):
        value = self._lookup.get(key)
        if value is None:
            value = self._method(key)
            self._lookup[key] = value
        return value
