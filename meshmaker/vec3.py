from .geometry import near, isnear, orient2d
from .quat import quat
import numpy as np


class vec3:

    @classmethod
    def O(cls):
        return cls(0, 0, 0)

    @classmethod
    def U(cls):
        return cls(1, 1, 1)

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

    def xy(self):
        return vec3(self.x, self.y, 0)

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

    def scl(self, o):
        if isinstance(o, vec3):
            return self.set(self.x * o.x, self.y * o.y, self.z * o.z)
        else:
            return self.set(self.x * o, self.y * o, self.z * o)

    def __repr__(self):
        return 'vec3({:.4f}, {:.4f}, {:.4f})'.format(self.x, self.y, self.z)

    def __init__(self, x, y, z):
        self.set(x, y, z)

    def set(self, x, y, z):
        self.x, self.y, self.z = x, y, z
        return self

    def near(self, o, e=0.00001):
        """Snap to other if sufficiently near."""
        raise NotImplementedError
        # TODO: this will snap for any axis despite the others...
        self.set(near(self.x, o.x, e),
                 near(self.y, o.y, e),
                 near(self.z, o.z, e))

    def isnear(self, o, e=0.00001):
        if isnear(self.x, o.x, e):
            if isnear(self.y, o.y, e):
                if isnear(self.z, o.z, e):
                    return True
        return False

    def axy(self, o):
        cosa = (self.dot(o) / (self.mag() * o.mag()))
        cosa = max(min(cosa, -1), 1)
        return np.arccos(cosa)

    def saxy(self, o):
        a = (np.arctan2(self.x, self.y) - np.arctan2(o.x, o.y))
        return a + 2 * np.pi if a < 0 else a

    def d(self, o):
        return (o - self).mag()

    def dxy(self, o):
        dx = o.x - self.x
        dy = o.y - self.y
        return np.sqrt(dx ** 2 + dy ** 2)

    def dexy(self, u, v):
        uv = (v - u)
        c = self.dot(uv)
        if c < u.dot(uv):
            return self.dxy(u)
        elif c > v.dot(uv):
            return self.dxy(v)
        else:
            nm = uv.crs(vec3.Z())
            return abs(self.dot(nm) - u.dot(nm))

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
        return self.set(x, y, z)

    def trn(self, o):
        return self.set(self.x + o.x, self.y + o.y, self.z + o.z)

    def crs(self, o):
        return vec3(self.y * o.z - self.z * o.y,
                    self.z * o.x - self.x * o.z,
                    self.x * o.y - self.y * o.x)

    def dot(self, o):
        return self.x * o.x + self.y * o.y + self.z * o.z

    def mag(self):
        return np.sqrt(self.dot(self))

    def nrm(self):
        mag = self.mag()
        if mag > 0:
            return vec3(self.x / mag, self.y / mag, self.z / mag)
        else:
            return vec3.O()

    def nrmd(self):
        mag = self.mag()
        if mag > 0:
            return self.set(self.x / mag, self.y / mag, self.z / mag)
        else:
            return self

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

    def onsxy(self, u, v, ie=False):
        #e = gtl.epsilon_c
        e = 0.00001
        if orient2d(self, u, v) == 0:
            tn = (v - u)
            a, b, c = self.dot(tn), u.dot(tn), v.dot(tn)
            try:
                b, c = (c, b) if c < b else (b, c)
            except:
                print(a, b, c)
                raise
            a = near(a, b, e)
            a = near(a, b, e)
            if b <= a and a <= c:
                if (a - b < e) or (c - a < e):
                    return (True if ie else False)

                #if (s1.x != s2.x): # S is not  vertical
                #if not gtl.isnear_c(s1.x,s2.x,gtl.epsilon_c): # S is not  vertical
                #    if (s1.x <= self.x and self.x <= s2.x):return 1
                #    if (s1.x >= self.x and self.x >= s2.x):return 1
                #else: # S is vertical, so test y  coordinate
                #    if (s1.y <= self.y and self.y <= s2.y):return 1
                #    if (s1.y >= self.y and self.y >= s2.y):return 1

        return False

    def inbxy(self, loop, ie=False):
        wn = 0
        for i in range(len(loop)):
            u, v = loop[i - 1], loop[i]
            if self.onsxy(u, v, ie=True):
                return 1 if ie else 0
            isleft = ((u.x - self.x) * (v.y - self.y) - (v.x - self.x) * (u.y - self.y))
            if u.y <= self.y:
                if v.y > self.y:
                    if isleft > 0:
                        wn += 1
            else:
                if v.y <= self.y:
                    if isleft < 0:
                        wn -= 1
        return wn

    def onbxy(self, loop, ie=False):
        for i in range(len(loop)):
            u, v = loop[i - 1], loop[i]
            if self.onsxy(u, v, ie=ie):
                return True
        return False

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


