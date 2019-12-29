from .base import Base
from .vec3 import vec3
from .quat import quat
import numpy as np


class TForm(Base):

    def __repr__(self):
        return f'tform:\n {self.position}\n {self.rotation}\n {self.scale}'

    def __init__(self,
                 position=None, rotation=None, scale=None,
                 children=None, parent=None, **kws):
        super().__init__(**kws)
        self.position = vec3.O() if position is None else position
        self.rotation = quat.O() if rotation is None else rotation
        self.scale = vec3.U() if scale is None else scale
        self.children = [] if children is None else children
        self.parent = parent
        if parent is not None:
            parent.children.append(self)

    def add(self, child):
        child.parent = self
        self.children.append(child)

    def _transform(self, world):
        local = (world * self.scale).rot(self.rotation) + self.position
        if self.parent is None:
            return local
        else:
            return self.parent._transform(local)

    def _inverse(self, local):
        world = (local + self.position.fp()).rot(self.rotation.fp()) * self.scale.inv()
        if self.parent is None:
            return world
        else:
            return self.parent._inverse(world)

    def transform(self, world):
        if isinstance(world, (list, tuple)):
            return list(self.transform(p) for p in world)
        else:
            return self._transform(world)

    def inverse(self, local):
        if isinstance(local, (list, tuple)):
            return list(self.inverse(p) for p in local)
        else:
            return self._inverse(local)

    def apply(self, f):
        f(self)
        for child in self.children:
            child.apply(f)
