import numpy as np
from meshmaker.vec3 import vec3
from meshmaker.quat import quat
from meshmaker.tform import TForm
from meshmaker.img import perlin, proximal, normalize
from meshmaker.geometry import slide, isnear


class field:
    """"""

    def __new__(cls, methods):
        def f(p):
            q = cls.zero()
            for m in methods:
                q = q + m(p)
            return q
        return f


    @staticmethod
    def decay(decay, weight):
        return lambda o, r: o * (weight * np.exp(-decay * r ** 2))


    @staticmethod
    def zero():
        raise NotImplementedError


class scalar_field(field):
    """"""

    @classmethod
    def zero(cls):
        return 0.0


    @classmethod
    def radial(cls, q, z, decay=0.0, weight=1.0):
        decay = cls.decay(decay, weight)
        def eigen(p):
            return decay(z, p.dxy(q))
        return eigen


    @classmethod
    def edge(cls, u, v, z, decay, weight):
        decay = cls.decay(decay, weight)
        def eigen(p):
            return decay(z, p.dexy(u, v))
        return eigen


    @classmethod
    def loop(cls, loop, z, decay, weight):
        loop = [cls.edge(u, v, z, decay, weight) for u, v in slide(loop, 2)]
        return cls(loop)


    @classmethod
    def topography(cls, img, tf, weight):
        def eigen(p):
            local = tf.transform(p)
            x = min(img.shape[1] - 1, int(local.x))
            y = min(img.shape[0] - 1, int(local.y))
            z = img[y, x]
            return z * weight
        return eigen


class image_field(scalar_field):
    """"""

    def __new__(cls, tf, img):
        es = (scalar_field.topography(img, tf, 1.0), )
        fd = super().__new__(cls, es)
        fd.img = img
        fd.tf = tf
        return fd


class height_field(image_field):
    """"""

    def __new__(cls, origin, radius, resolution, landmasses,
                dz=1, noise=False, prox_f=None):
        pixels = int(2 * radius * resolution)
        tf = cls.wtoi(origin, radius, pixels)
        n = perlin(pixels, 16 * 8, 8, 4) if noise else 1.0
        f = (lambda d: max(0, d) ** 0.8) if prox_f is None else prox_f
        p = [proximal(pixels, tf.transform(lm), f) for lm in landmasses]
        p = sum(p)
        p = normalize(p)
        height = n * p * dz

		#import cv2
        #for j in range(2):
        #    height = cv2.blur(height, (8, 8))

        return super().__new__(cls, tf, height)

    @staticmethod
    def wtoi(origin, radius, resolution):
        dp = origin.tov(vec3(resolution / 2, resolution / 2, 0))
        ds = vec3(resolution / (2 * radius), resolution / (2 * radius), 1)
        tf = TForm(dp, None, ds)
        return tf


class vec3_field(field):
    """"""

    @staticmethod
    def zero():
        return vec3.O()


    @classmethod
    def radial(cls, q, decay, weight):
        decay = cls.decay(decay, weight)
        def eigen(p):
            return decay((p - q).nrm(), p.dxy(q))
        return eigen


    @classmethod
    def grid(cls, q, major, decay, weight):
        decay = cls.decay(decay, weight)
        def eigen(p):
            return decay(major, p.dxy(q))
        return eigen


    @classmethod
    def edge(cls, u, v, decay, weight):
        decay = cls.decay(decay, weight)
        tangent = (v - u).nrm()
        #normal = tangent.crs(vec3.Z()).nrm()
        normal = tangent.crs(vec3.Z()).nrm().fp()
        def eigen(p):
            isleft = ((u.x - p.x) * (v.y - p.y) - (v.x - p.x) * (u.y - p.y))
            return decay(normal * (1.0 if isleft else -1.0), p.dexy(u, v))
            #return decay(normal * (-1 if isleft else 1), p.dexy(u, v))
            #return decay(tangent, p.dexy(u, v)) 
        return eigen


    @classmethod
    def loop(cls, loop, decay, weight):
        loop = [cls.edge(u, v, decay, weight) for u, v in slide(loop, 2)]
        return cls(loop)


    @classmethod
    def topography(cls, terrain, weight):
        dx, dy = np.gradient(terrain.img)
        dx = image_field(terrain.tf, dx)
        dy = image_field(terrain.tf, dy)
        def eigen(p):
            dpdx = dx(p)
            dpdy = dy(p)
            return vec3(dpdx, dpdy, 0) * weight * max(1, np.sqrt(dpdx ** 2 + dpdy ** 2))
        return eigen


class trace_field(vec3_field):
    """"""

    def __new__(cls, methods, ds=1, alpha=0):
        fd = super().__new__(cls, methods)
        def f(p):
            """trace a sequence of new edges using fields"""
            last_fd = fd(p).nrm()

            while True:

                next_fd = fd(p).nrm()
                sa = last_fd.saxy(next_fd)
                last_fd = next_fd
                if isnear(sa, np.pi):
                    print('singularity flip')
                    next_fd = next_fd.fp()

                q = p + (next_fd * ds).rot(quat.av(alpha, vec3.Z()))
                yield q
                p = q
        return f
