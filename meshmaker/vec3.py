from .geometry import near, isnear, orient2d
from .quat import quat
import numpy as np


class vec3:

    @classmethod
    def O(cls):
        return cls(0, 0, 0)

    @classmethod
    def X(cls):
        return cls(1, 0, 0)

    @classmethod
    def Y(cls):
        return cls(0, 1, 0)

    @classmethod
    def Z(cls):
        return cls(0, 0, 1)

    @classmethod
    def com(cls, pts):
        xs, ys, zs = zip(*[(p.x, p.y, p.z) for p in pts])
        return cls(sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))

    def cp(self):
        return vec3(self.x, self.y, self.z)

    def fp(self):
        return vec3(-self.x, -self.y, -self.z)

    def inv(self):
        return vec3(1 / self.x, 1 / self.y, 1 / self.z)

    def __add__(self, o):
        if isinstance(o, vec3):
            return vec3(self.x + o.x, self.y + o.y, self.z + o.z)
        else:
            return vec3(self.x + o, self.y + o, self.z + o)

    def __sub__(self, o):
        if isinstance(o, vec3):
            return vec3(self.x - o.x, self.y - o.y, self.z - o.z)
        else:
            return vec3(self.x - o, self.y - o, self.z - o)

    def __mul__(self, o):
        if isinstance(o, vec3):
            return vec3(self.x * o.x, self.y * o.y, self.z * o.z)
        else:
            return vec3(self.x * o, self.y * o, self.z * o)

    def __repr__(self):
        return 'vec3({:.4f}, {:.4f}, {:.4f})'.format(self.x, self.y, self.z)

    def __init__(self, x, y, z):
        self.set(x, y, z)

    def set(self, x, y, z):
        self.x, self.y, self.z = x, y, z

    def near(self, o, e=0.00001):
        """Snap to other if sufficiently near."""
        raise NotImplementedError
        self.set(near(self.x, o.x, e),
                 near(self.y, o.y, e),
                 near(self.z, o.z, e))

    def isnear(self, o, e=0.00001):
        x = isnear(self.x, o.x, e)
        y = isnear(self.y, o.y, e)
        z = isnear(self.z, o.z, e)
        return x and y and z

    def saxy(self, o):
        a = (np.arctan2(self.x, self.y) - np.arctan2(o.x, o.y))
        return a + 2 * np.pi if a < 0 else a

    def tov(self, o):
        return o - self

    def rot(self, q):
        x = self.dot(vec3(
            q.w ** 2 + q.x ** 2 - q.y ** 2 - q.z ** 2,
            2 * (q.x * q.y - q.w * q.z),
            2 * (q.x * q.z + q.w * q.y)))
        y = self.dot(vec3(
            2 * (q.x * q.y + q.w * q.z),
            q.w ** 2 - q.x ** 2 + q.y ** 2 - q.z ** 2,
            2 * (q.y * q.z - q.w * q.x)))
        z = self.dot(vec3(
            2 * (q.x * q.z - q.w * q.y),
            2 * (q.y * q.z + q.w * q.x),
            q.w ** 2 - q.x ** 2 - q.y ** 2 + q.z ** 2))
        self.x = x
        self.y = y
        self.z = z
        return self

    def trn(self, o):
        self.x += o.x
        self.y += o.y
        self.z += o.z
        return self

    def crs(self, o):
        return vec3(self.y * o.z - self.z * o.y,
                    self.z * o.x - self.x * o.z,
                    self.x * o.y - self.y * o.x)

    def dot(self, o):
        return self.x * o.x + self.y * o.y + self.z * o.z

    def mag(self):
        return np.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def nrm(self):
        mag = self.mag()
        return vec3(self.x / mag, self.y / mag, self.z / mag)

    def lerp(self, o, d):
        return self + (self.tov(o) * d)

    def ring(self, r, n, inscribe=True):
        alpha = np.pi * (2.0 / n)
        sr = r if inscribe else r / np.cos(alpha / 2.0)
        z = vec3(0, 0, 1)
        loop = []
        for x in range(n):
            loop.append(vec3(sr, 0, 0).rot(quat.av(x * alpha - alpha / 2.0, z)))
        return [p.trn(self) for p in loop]

    def line(self, o, n):
        line = []
        for i in range(n):
            t = (i + 1) / (n + 1)
            line.append(self.lerp(o, t))
        return line

    def fan(self, r, n, inscribe=True):
        ring = self.ring(r, n, inscribe)
        return [(self, ring[i - 1], ring[i]) for i in range(n)]

    def insxy(self, u, v):
        if isnear(orient2d(self, u, v), 0):
            uv = u.tov(v)
            u_uv = u.dot(uv)
            v_uv = v.dot(uv)
            self_uv = self.dot(uv)
            if v_uv < u_uv:
                u_uv, v_uv = v_uv, u_uv
            return u_uv <= self_uv and self_uv <= v_uv
        else:
            return False

        """
        e = 0.01
        if gtl.orient2d_c(self, s1, s2) == 0:
            tn = s1.tov(s2)
            selfprj, s1prj, s2prj = self.dot(tn), s1.dot(tn), s2.dot(tn)
            if s2prj < s1prj:
                s1prj, s2prj = s2prj, s1prj
            selfprj = gtl.near_c(selfprj, s1prj, e)
            selfprj = gtl.near_c(selfprj, s2prj, e)
            if s1prj <= selfprj and selfprj <= s2prj:
                #if self.isnear(s1) or self.isnear(s2):
                if (selfprj - s1prj < e) or (s2prj - selfprj < e):
                    return (1 if ie else 0)

                #if (s1.x != s2.x): # S is not  vertical
                #if not gtl.isnear_c(s1.x,s2.x,gtl.epsilon_c): # S is not  vertical
                #    if (s1.x <= self.x and self.x <= s2.x):return 1
                #    if (s1.x >= self.x and self.x >= s2.x):return 1
                #else: # S is vertical, so test y  coordinate
                #    if (s1.y <= self.y and self.y <= s2.y):return 1
                #    if (s1.y >= self.y and self.y >= s2.y):return 1

        return 0
        """
        raise NotImplementedError

    def onsxy(self, u, v):
        raise NotImplementedError

    def inbxy(self, loop):
        """
        cdef int wn = 0
        cdef int px
        cdef int pcnt = len(ps)
        cdef vec3 p1,p2
        cdef float isleft
        for px in range(pcnt):
            p1,p2 = ps[px-1],ps[px]
            if self.onsxy_c(p1,p2,1):return 0
            isleft = ((p1.x-self.x)*(p2.y-self.y)-(p2.x-self.x)*(p1.y-self.y))
            if p1.y <= self.y:
                if p2.y > self.y:
                    if isleft > 0:
                        wn += 1
            else:
                if p2.y <= self.y:
                    if isleft < 0:
                        wn -= 1
        return wn
        """
        raise NotImplementedError

    def onbxy(self, loop):
        raise NotImplementedError

    def incircle(self, a, b, c):
        """True if self is inside the circumcircle of the triangle a, b, c"""
        m11, m12 = a.x - d.x, a.y - d.y
        m13 = m11 * m11 + m12 * m12
        m21, m22 = b.x - d.x, b.y - d.y
        m23 = m21 * m21 + m22 * m22
        m31, m32 = c.x - d.x, c.y - d.y
        m33 = m31 * m31 + m32 * m32
        det1 = m11 * (m22 * m33 - m23 * m32)
        det2 = m12 * (m21 * m33 - m23 * m31)
        det3 = m13 * (m21 * m32 - m22 * m31)
        return near(det1 - det2 + det3, 0)

    def __orient2d(self, a, b):
        """Signed area of the triangle a - self, b - self
            0 if a, b, and self are colinear
        """
        m11 = a.x - self.x
        m12 = a.y - self.y
        m21 = b.x - self.x
        m22 = b.y - self.y
        det = m11 * m22 - m12 * m21
        return near(det, 0)


if __name__ == '__main__':
    def show(u, v):
        a = u.saxy(v)
        print(u, v, (180 / np.pi) * a)
    show(vec3(-1, 0, 0), vec3(1, 1, 0))
    show(vec3(1, 1, 0), vec3(-1, 0, 0))
    show(vec3( 1, 0, 0), vec3(1, 1, 0))
    show(vec3( 1, 0, 0), vec3(1, 0, 0))
    show(vec3( 1, -1, 0), vec3(-1, 1, 0))


