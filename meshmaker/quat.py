from .vec3 import vec3
from .geometry import isnear
import numpy as np


class quat:

    def __repr__(self):
        return f'quat({self.w:.4f}, {self.x:.4f}, {self.y:.4f}, {self.z:.4f})'

    def __init__(self, w, x, y, z):
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def O(cls):
        return cls(1, 0, 0, 0)

    @classmethod
    def av(cls, a, v):
        v = v * (np.sin(a / 2.0) / v.mag())
        return cls(np.cos(a / 2.0), v.x, v.y, v.z)

    @classmethod
    def uu(cls, x, y):
        a = x.ang(y)
        if a == 0.0:
            return cls(1, 0, 0, 0)
        else:
            v = x.crs(y).nrm()
            if isnear(v.dot(v), 0):
                v = x.crs(vec3.X()).nrm()
                if isnear(v.dot(v), 0):
                    v = x.crs(vec3.Y()).nrm()
                    if isnear(v.dot(v), 0):
                        raise
                #v = vec3.Z()
            return cls.av(a, v)

    @classmethod
    def toxy(cls, v):
        vz = v.nrm().z
        if isnear(vz, -1):
            return cls(0, 1, 0, 0)
        elif not isnear(vz, 1):
            return cls.uu(v, vec3.Z())
        else:
            return cls.av(0, vec3.Z())

    @classmethod
    def rotz(cls, a):
        return cls.av(a, vec3.Z())

    def fp(self):
        return quat(-self.w, self.x, self.y, self.z)

    def rot(self, ps):
        for p in ps:
            p.rot(self)
        return ps


