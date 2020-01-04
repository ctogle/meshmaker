from .vec3 import vec3
from .quat import quat
import numpy as np


class mat44:

    @classmethod
    def I(cls):
        return cls(1, 0, 0, 0,
                   0, 1, 0, 0,
                   0, 0, 1, 0,
                   0, 0, 0, 1)

    @classmethod
    def tform(cls, t=None, q=None, s=None):
        if t is None:
            m14 = m24 = m34 = 0
        else:
            m14, m24, m34 = t.x, t.y, t.z
        if q is None:
            m11 = m22 = m33 = 1
            m12 = m13 = m23 = m21 = m31 = m32 = 0
        else:
            v = vec3(q.x, q.y, q.z).nrm()
            a = 2 * np.arccos(q.w)
            cosa, sina = np.cos(a), np.sin(a)
            m11 = v.x * v.x * (1 - cosa) + cosa
            m12 = v.y * v.x * (1 - cosa) - v.z * sina
            m13 = v.z * v.x * (1 - cosa) + v.y * sina
            m21 = v.x * v.y * (1 - cosa) + v.z * sina
            m22 = v.y * v.y * (1 - cosa) + cosa
            m23 = v.z * v.y * (1 - cosa) - v.x * sina
            m31 = v.x * v.z * (1 - cosa) - v.y * sina
            m32 = v.y * v.z * (1 - cosa) + v.x * sina
            m33 = v.z * v.z * (1 - cosa) + cosa
        if s is not None:
            m11 *= s.x
            m22 *= s.y
            m33 *= s.z
        return cls(m11, m12, m13, m14,
                   m21, m22, m23, m24,
                   m31, m32, m33, m34,
                     0,   0,   0,   1)

    def __repr__(self):
        return f'{self.m11} {self.m12} {self.m13} {self.m14}\
               \n{self.m21} {self.m22} {self.m23} {self.m24}\
               \n{self.m31} {self.m32} {self.m33} {self.m34}\
               \n{self.m41} {self.m42} {self.m43} {self.m44}'

    def __init__(self,
                 m11, m12, m13, m14,
                 m21, m22, m23, m24,
                 m31, m32, m33, m34,
                 m41, m42, m43, m44):
        self.m11, self.m12, self.m13, self.m14 = m11, m12, m13, m14
        self.m21, self.m22, self.m23, self.m24 = m21, m22, m23, m24
        self.m31, self.m32, self.m33, self.m34 = m31, m32, m33, m34
        self.m41, self.m42, self.m43, self.m44 = m41, m42, m43, m44

    def __mul__(self, o):
        if isinstance(o, vec3):
            return vec3(
                self.m11 * o.x + self.m12 * o.y + self.m13 * o.z + self.m14,
                self.m21 * o.x + self.m22 * o.y + self.m23 * o.z + self.m24,
                self.m31 * o.x + self.m32 * o.y + self.m33 * o.z + self.m34)
        elif isinstance(o, mat44):
            return mat44(
                self.m11 * o.m11 + self.m12 * o.m21 + self.m13 * o.m31 + self.m14 * o.m41,
                self.m11 * o.m12 + self.m12 * o.m22 + self.m13 * o.m32 + self.m14 * o.m42,
                self.m11 * o.m13 + self.m12 * o.m23 + self.m13 * o.m33 + self.m14 * o.m43,
                self.m11 * o.m14 + self.m12 * o.m24 + self.m13 * o.m34 + self.m14 * o.m44,
                self.m21 * o.m11 + self.m22 * o.m21 + self.m23 * o.m31 + self.m24 * o.m41,
                self.m21 * o.m12 + self.m22 * o.m22 + self.m23 * o.m32 + self.m24 * o.m42,
                self.m21 * o.m13 + self.m22 * o.m23 + self.m23 * o.m33 + self.m24 * o.m43,
                self.m21 * o.m14 + self.m22 * o.m24 + self.m23 * o.m34 + self.m24 * o.m44,
                self.m31 * o.m11 + self.m32 * o.m21 + self.m33 * o.m31 + self.m34 * o.m41,
                self.m31 * o.m12 + self.m32 * o.m22 + self.m33 * o.m32 + self.m34 * o.m42,
                self.m31 * o.m13 + self.m32 * o.m23 + self.m33 * o.m33 + self.m34 * o.m43,
                self.m31 * o.m14 + self.m32 * o.m24 + self.m33 * o.m34 + self.m34 * o.m44,
                self.m41 * o.m11 + self.m42 * o.m21 + self.m43 * o.m31 + self.m44 * o.m41,
                self.m41 * o.m12 + self.m42 * o.m22 + self.m43 * o.m32 + self.m44 * o.m42,
                self.m41 * o.m13 + self.m42 * o.m23 + self.m43 * o.m33 + self.m44 * o.m43,
                self.m41 * o.m14 + self.m42 * o.m24 + self.m43 * o.m34 + self.m44 * o.m44)
        else:
            print(type(o), o)
            raise NotImplementedError
