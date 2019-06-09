""""""
import numpy as np
from .vec3 import vec3


class scalar_field:
    """"""

    @staticmethod
    def apply(img, tf, world):
        local = tf.transform(world)
        x = min(img.shape[1] - 1, int(local.x))
        y = min(img.shape[0] - 1, int(local.y))
        z = img[y, x]
        return z

    def __init__(self, tf, img):
        self.tf = tf
        self.img = img

    def __call__(self, p):
        return self.apply(self.img, self.tf, p)


class vec3_field:
    """"""

    @staticmethod
    def decay(decay, weight):
        return lambda o, r: o * (weight * np.exp(-decay * r ** 2))

    @classmethod
    def radial(cls, q, decay, weight):
        decay = cls.decay(decay, weight)
        def eigen(p):
            return decay((p - q).nrm(), p.dxy(q))
        return eigen

    @classmethod
    def edge(cls, u, v, decay, weight):
        decay = cls.decay(decay, weight)
        tangent = (v - u).nrm()
        normal = tangent.crs(vec3.Z()).nrm()
        def eigen(p):
            isleft = ((u.x - p.x) * (v.y - p.y) - (v.x - p.x) * (u.y - p.y))
            return decay(normal * (1 if isleft else -1), p.dexy(u, v))
            #return decay(tangent, p.dexy(u, v))
        return eigen

    @classmethod
    def grid(cls, q, major, decay, weight):
        decay = cls.decay(decay, weight)
        def eigen(p):
            return decay(major, p.dxy(q))
        return eigen

    @classmethod
    def topography(cls, terrain, weight):
        dx, dy = np.gradient(terrain.img)
        dx = scalar_field(terrain.tf, dx)
        dy = scalar_field(terrain.tf, dy)
        def eigen(p):
            dpdx = dx(p)
            dpdy = dy(p)
            return vec3(dpdx, dpdy, 0) * weight * max(1, np.sqrt(dpdx ** 2 + dpdy ** 2))
        return eigen

    @classmethod
    def parse_elements(cls, element, decay, weight):
        if isinstance(element, vec3):
            yield cls.radial(element, decay, weight)
        elif isinstance(element, tuple):
            yield cls.grid(element[0], element[1], decay, weight)
        elif isinstance(element, list):
            for u, v in slide(element, 2):
                yield cls.edge(u, v, decay, weight)
        elif isinstance(element, scalar_field):
            yield cls.topography(element, weight)
        else:
            raise ValueError

    @staticmethod
    def apply(world, elements):
        v = vec3.O()
        for e in elements:
            v += e(world)
        return v

    def __new__(cls, elements):
        methods = []
        for e in elements:
            methods.extend(list(cls.parse_elements(*e)))
        return lambda p: cls.apply(p, methods)


