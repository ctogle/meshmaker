from .geometry import near, isnear, orient2d, inrng, slide
import numpy as np
import math


class vec3:

    @classmethod
    def O(cls):
        return cls(0, 0, 0)

    @classmethod
    def U(cls, r=1):
        return cls(r, r, r)

    @classmethod
    def X(cls, r=1):
        return cls(r, 0, 0)

    @classmethod
    def Y(cls, r=1):
        return cls(0, r, 0)

    @classmethod
    def Z(cls, r=1):
        return cls(0, 0, r)

    @classmethod
    def nX(cls, r=1):
        return cls(-r, 0, 0)

    @classmethod
    def nY(cls, r=1):
        return cls(0, -r, 0)

    @classmethod
    def nZ(cls, r=1):
        return cls(0, 0, -r)

    @classmethod
    def com(cls, pts):
        xs, ys, zs = zip(*[(p.x, p.y, p.z) for p in pts])
        return cls(sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))

    def isnan(self):
        return math.isnan(self.x) or math.isnan(self.y) or math.isnan(self.z)

    def isO(self):
        return self.x == 0 and self.y == 0 and self.z == 0

    def xy(self):
        return vec3(self.x, self.y, 0)

    def yz(self):
        return vec3(0, self.y, self.z)

    def zx(self):
        return vec3(self.x, 0, self.z)

    def cp(self):
        return vec3(self.x, self.y, self.z)

    def fp(self):
        return vec3(-self.x, -self.y, -self.z)

    def inv(self):
        return vec3(1 / self.x, 1 / self.y, 1 / self.z)

    @classmethod
    def sum(cls, pts):
        r = vec3.O()
        for p in pts:
            r.trn(p)
        return r

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

    def __iter__(self):
        yield self.x
        yield self.y
        yield self.z

    def scl(self, o):
        if isinstance(o, vec3):
            return self.set(self.x * o.x, self.y * o.y, self.z * o.z)
        else:
            return self.set(self.x * o, self.y * o, self.z * o)
    def sclps(self, os):
        for o in os:
            o.set(self.x * o.x, self.y * o.y, self.z * o.z)
        return os

    def quant(self, n=4):
        return vec3(round(self.x, 4), round(self.y, 4), round(self.z, 4))

    def __repr__(self):
        return 'vec3({:.4f}, {:.4f}, {:.4f})'.format(self.x, self.y, self.z)

    def __init__(self, x, y, z):
        self.set(x, y, z)

    def setto(self, o):
        return self.set(o.x, o.y, o.z)

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

    def ang(self, o):
        cosa = (self.dot(o) / (self.mag() * o.mag()))
        cosa = min(max(cosa, -1), 1)
        return np.arccos(cosa)

    def axy(self, o):
        return self.xy().ang(o.xy())

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

    def dlxy(self, loop):
        return min([self.dexy(u, v) for u, v in slide(loop, 2)])

    def tov(self, o):
        return o - self

    def rot(self, q):
        x = self.x * (q.w ** 2 + q.x ** 2 - q.y ** 2 - q.z ** 2) +\
            self.y * (2 * (q.x * q.y - q.w * q.z)) +\
            self.z * (2 * (q.x * q.z + q.w * q.y))
        y = self.x * (2 * (q.x * q.y + q.w * q.z)) +\
            self.y * (q.w ** 2 - q.x ** 2 + q.y ** 2 - q.z ** 2) +\
            self.z * (2 * (q.y * q.z - q.w * q.x))
        z = self.x * (2 * (q.x * q.z - q.w * q.y)) +\
            self.y * (2 * (q.y * q.z + q.w * q.x)) +\
            self.z * (q.w ** 2 - q.x ** 2 - q.y ** 2 + q.z ** 2)
        return self.set(x, y, z)

    def trn(self, o):
        return self.set(self.x + o.x, self.y + o.y, self.z + o.z)

    def xtrn(self, dx):
        return self.set(self.x + dx, self.y, self.z)

    def ytrn(self, dy):
        return self.set(self.x, self.y + dy, self.z)

    def ztrn(self, dz):
        return self.set(self.x, self.y, self.z + dz)

    def trnps(self, os, cp=False):
        for o in os:
            o.set(self.x + o.x, self.y + o.y, self.z + o.z)
        return os

    def trnpy(self, py):
        return (self.trnps(py[0]), [self.trnps(h) for h in py[1]])

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
        from .quat import quat
        alpha = np.pi * (2.0 / n)
        sr = r if inscribe else r / np.cos(alpha / 2.0)
        z = vec3(0, 0, 1)
        loop = []
        for x in range(n):
            loop.append(vec3(sr, 0, 0).rot(quat.av(x * alpha - alpha / 2.0, z)))
        return [p.trn(self) for p in loop]

    def line(self, o, n, ends=False):
        line = []
        if ends:
            line.append(self.cp())
        for i in range(n):
            t = (i + 1) / (n + 1)
            line.append(self.lerp(o, t))
        if ends:
            line.append(o.cp())
        return line

    def spline(self, o, st, ot, n, alpha=0.5):
        #n += 1
        n -= 1
        ps = [self, self + st, o + ot, o]
        x, y, z = zip(*[(p.x, p.y, p.z) for p in ps])
        t = np.cumsum([0] + [u.d(v) ** alpha for u, v in slide(ps, 2, 1)])
        x = self.catmull(x, t, n)[1:-1]
        y = self.catmull(y, t, n)[1:-1]
        z = self.catmull(z, t, n)[1:-1]
        return [vec3(x, y, z) for x, y, z in zip(x, y, z)]

    def catmull(self, xs, t, n):

        def coordinate(xs, t, k):
            l01 = xs[0]*(t[1] - k)/(t[1] - t[0]) + xs[1]*(k - t[0])/(t[1] - t[0])
            l12 = xs[1]*(t[2] - k)/(t[2] - t[1]) + xs[2]*(k - t[1])/(t[2] - t[1])
            l23 = xs[2]*(t[3] - k)/(t[3] - t[2]) + xs[3]*(k - t[2])/(t[3] - t[2])
            l012 = l01*(t[2] - k)/(t[2] - t[0]) + l12*(k - t[0])/(t[2] - t[0])
            l123 = l12*(t[3] - k)/(t[3] - t[1]) + l23*(k - t[1])/(t[3] - t[1])
            c12 = l012*(t[2] - k)/(t[2] - t[1]) + l123*(k - t[1])/(t[2] - t[1])
            return c12

        curve = [xs[0]]
        for i in range(1, len(xs) - 2):
            for j in range(n):
                k = t[1] + (j / n) * (t[2] - t[1])
                x = coordinate(xs[i - 1:i + 3], t[i - 1:i + 3], k)
                curve.append(x)
        curve.append(xs[-2])
        curve.append(xs[-1])
        return curve

    def fan(self, r, n, inscribe=True):
        ring = self.ring(r, n, inscribe)
        return [(self, ring[i - 1], ring[i]) for i in range(n)]

    def insxy(self, u, v, e=0.00001):
        return self.onsxy(u, v, ie=False, e=e)

        '''
        if isnear(orient2d(self, u, v, epsilon=e), 0):
            uv = u.tov(v)
            u_uv = u.dot(uv)
            v_uv = v.dot(uv)
            self_uv = self.dot(uv)
            if v_uv < u_uv:
                u_uv, v_uv = v_uv, u_uv
            return u_uv <= self_uv and self_uv <= v_uv
        else:
            return False
        '''

    def onsxy(self, u, v, ie=False, e=0.00001):
        perp = self.leftof(u, v, e)
        if perp == 0:
            t = (v - u).nrm()
            a, b, c = u.dot(t), self.dot(t), v.dot(t)
            return (ie and (a == b or b == c)) or inrng(b, a, c)
        else:
            return False

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

    def leftof(self, u, v, e=0.00001):
        """compute dxy to line containing edge uv"""
        n = vec3.Z().crs(v - u).nrm()
        return near(self.dot(n) - u.dot(n), 0, e)

    def inbxy(self, loop, ie=False, e=0.00001):
        wn = 0
        for i in range(len(loop)):
            u, v = loop[i - 1], loop[i]
            if self.onsxy(u, v, ie=True, e=e):
                return 1 if ie else 0
            x = (u.x - self.x) * (v.y - self.y)
            y = (v.x - self.x) * (u.y - self.y)
            isleft = near(x - y, 0)
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
        m11, m12 = a.x - self.x, a.y - self.y
        m13 = m11 * m11 + m12 * m12
        m21, m22 = b.x - self.x, b.y - self.y
        m23 = m21 * m21 + m22 * m22
        m31, m32 = c.x - self.x, c.y - self.y
        m33 = m31 * m31 + m32 * m32
        det1 = m11 * (m22 * m33 - m23 * m32)
        det2 = m12 * (m21 * m33 - m23 * m31)
        det3 = m13 * (m21 * m32 - m22 * m31)
        return near(det1 - det2 + det3, 0)

