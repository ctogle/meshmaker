from .base import Base
from .vec3 import vec3
from .quat import quat
from .mat44 import mat44
import numpy as np


class TForm(Base):

    def world(self):
        if self.parent is None:
            base = mat44.I()
        else:
            base = self.parent.world()
        return base * self.transformation

    def transform(self, vs):
        w = self.world()
        if isinstance(vs, (list, tuple)):
            return [w * v for v in vs]
        else:
            return w * vs

    def __init__(self, t=None, q=None, s=None,
                 parent=None, children=None, **kws):
        super().__init__(**kws)
        self.transformation = mat44.tform(t, q, s)

        self.parent = parent
        if parent:
            parent.add(self)
        self.children = []
        if children is not None:
            for child in children:
                self.add(child)

    def add(self, child):
        child.parent = self
        if not child in self.children:
            self.children.append(child)




class ___TForm(Base):

    def __repr__(self):
        return f'tform:\n {self.position}\n {self.rotation}\n {self.scale}'

    def __init__(self,
                 position=None, rotation=None, scale=None,
                 children=None, parent=None, **kws):
        super().__init__(**kws)
        self.position = position
        self.rotation = rotation
        self.scale = scale
        #self.position = vec3.O() if position is None else position
        #self.rotation = quat.O() if rotation is None else rotation
        #self.scale = vec3.U() if scale is None else scale
        self.children = []
        if children is not None:
            for child in children:
                self.add(child)
        self.parent = parent
        if parent is not None:
            parent.add(self)

    def add(self, child):
        child.parent = self
        self.children.append(child)

    def _transform(self, world):
        local = world.cp()
        if self.scale:
            local.scl(self.scale)
        if self.rotation:
            local.rot(self.rotation)
        if self.position:
            local.trn(self.position)
        #local = (world * self.scale).rot(self.rotation) + self.position
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
