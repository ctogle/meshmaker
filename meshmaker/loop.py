from .vec3 import vec3
from .quat import quat
from .geometry import slide, isnear, sintsxyp, bbox, loop_area, loop_normal
from .planargraph import planargraph
from .delaunay import triangulation
from .plt import *
import numpy as np
from contextlib import contextmanager


class loops:

    @property
    def N(self):
        """Compute the normal vector of self"""
        N = []
        for ps in self.loops:
            N.append(loop_normal(ps))
            #N.append(vec3.O())
            #for u, v, w in slide(ps, 3):
            #    uv, vw = (v - u), (w - v)
            #    N[-1].trn(uv.crs(vw))
            #assert not N[-1].isO(), 'zero vector in loop.N ...'
        return vec3.sum(N).nrm()

    def _toxy(self, *loops):
        if self._q is None:
            N = self.N
            if N.dot(N) > 0:
                q = quat.toxy(N)
                if not isnear(q.w, 0):
                    self.rot(q)
                    for loop in loops:
                        q.rot(loop)
                self._q = q
        else:
            q = self._q
            if not isnear(q.w, 0):
                for loop in loops:
                    q.rot(loop)
        return self

    def _fromxy(self, *loops):
        if not self._q is None:
            q = self._q.fp()
            if not isnear(q.w, 0):
                self.rot(q)
                for loop in loops:
                    q.rot(loop)
            self._q = None
        return self

    @property
    def area(self):
        """Compute the area contained by self"""
        self._toxy()
        area = 0.0
        for ps in self.loops:
            area += loop_area(ps)
        for ps in self.holes:
            area -= loop_area(ps)
        self._fromxy()
        return area

    def contains(self, p):
        """Check if p is on the interior/boundary of self"""
        self._toxy()
        for loop in self.loops:
            if p.inbxy(loop, True):
                return True
        self._fromxy()
        return False

    def plot(self, ax, toxy=True, **kws):
        if toxy:
            self._toxy()
        for ps in self.loops:
            plot_loop(ax, ps, **kws)
        for ps in self.holes:
            kws['ls'] = kws.get('ls', '--')
            plot_loop(ax, ps, **kws)
        if toxy:
            self._fromxy()

    def triangulation(self, e=0.0001, h=None, r=10000):
        """Compute a delaunay triangulation of self"""
        # TODO: handle loop hierarchy...
        assert len(self.loops) == 1
        self._toxy()
        py = (self.loops[0], [h.ps for h in self.holes])
        t = triangulation(py, e, h, r)
        q.rot(t.points)
        self._fromxy()
        return t

    def __init__(self, loops=(), holes=()):
        self.loops = []
        self.holes = []
        for l in loops:
            self.al(l)
        for h in holes:
            self.ah(h)
        self._q = None

    def al(self, loop):
        self.loops.append([p.quant() for p in loop])

    def ah(self, hole):
        self.holes.append([p.quant() for p in hole])

    def _findparent(self, hole, mode='strict'):
        # must fit within exactly one parent contour
        # without intersecting that parent contour
        onb = not mode == 'strict'
        for c, l in enumerate(self.loops):
            #if all(p.inbxy(l, False) for p in hole):
            if any(p.inbxy(l, onb) for p in hole):
                parent = c
                break
        else:
            return
        for u, v in slide(self.loops[parent], 2):
            for x, y in slide(hole, 2):
                if sintsxyp(u, v, x, y):
                    return None if mode == 'strict' else parent
        return parent

    def _strict_embed(self, hole):
        # must not overlap/intersect any existent holes
        for h in self.holes:
            if any(p.inbxy(h, True) for p in hole):
                return
            for u, v in slide(h, 2):
                for x, y in slide(hole, 2):
                    if sintsxyp(u, v, x, y):
                        return
        parent = self._findparent(hole, mode='strict')
        # this hole can be safely added
        if parent is not None:
            self.ah(hole)
            return True
        return False

    def _clip_embed(self, hole):
        strictly_embedded = self._strict_embed(hole)
        if not strictly_embedded:
            parent = self._findparent(hole, mode='clip')
            if parent is not None:

                print('clip self to hole')
                #print(self.loops[parent])
                #print(hole)

                # check if it can be added non-strictly...
                pcp = self.__class__((self.loops[parent], ))
                hcp = self.__class__((hole, ))
                #qO = quat.O()
                #pcp._q = qO
                #hcp._q = qO
                #pcp._q = self._q
                #hcp._q = self._q
                reparent = pcp.difference(hcp, inplace=False)
                unparent = self.loops.pop(parent)
                # these new loops end up in the xy plane regardless of
                # orientation pre-op
                for newparent in reparent.loops:
                    self.loops.insert(parent, newparent)

                return True

        return False

    def embed(self, hole, mode='strict'):
        """Smarter self.ah"""
        hole = [p.quant() for p in hole]
        self._toxy(hole)
        if mode == 'strict':
            embedded = self._strict_embed(hole)
        elif mode == 'clip':
            embedded = self._clip_embed(hole)
        else:
            raise
        self._fromxy()
        #return self
        return embedded

    @classmethod
    def from_polygon(cls, loop, holes=None):
        #if holes is not None:
        #    holes = [cls(hole) for hole in holes]
        return cls([loop], holes)

    def rot(self, q):
        for ps in self.loops:
            q.rot(ps)
        for ps in self.holes:
            q.rot(ps)
        return self

    def cp(self):
        loops = [[p.cp() for p in ps] for ps in self.loops]
        holes = [[p.cp() for p in ps] for ps in self.holes]
        return self.__class__(loops, holes)

    def lexico(self, key=(lambda t: (t[1].z, t[1].x, t[1].y))):
        """Find the ordering of self.ps where the first two points
        have the lowest z-values possible"""
        # TODO: clean this up

        raise

        z = bbox(self.ps)[0].z
        c = 0
        for o, p in enumerate(self.ps):
            if isnear(p.z, z):
                if c:
                    break
                else:
                    c += 1
        o -= 1

        # o is the first p in self.ps with minimal z coord

        #o, _ = sorted(list(enumerate(self.ps)), key=key)[0]
        return [p.cp() for p in (self.ps[o:] + self.ps[:o])]

    def fp(self):
        for ps in self.loops:
            ps.reverse()
        for ps in self.holes:
            ps.reverse()
        return self

    def smooth(self, weight=0.5, iterations=1):
        """Smooth towards 1-ring neighbors in-place"""
        # TODO: handle holes
        for ps in self.loops:
            for j in range(iterations):
                dps = []
                for u, v, w in slide(ps, 3):
                    dp = (vec3.com((u, w)) - v) * weight
                    dps.append(dp)
                dps.insert(0, dps.pop(-1))
                for p, dp in zip(ps, dps):
                    p.trn(dp)
        return self

    def edgesplit(self, maxlength):
        """Compute new loop where no edge exceeds maxlength"""
        # TODO: handle holes
        loops = []
        for ps in self.loops:
            points = []
            for u, v in slide(ps, 2):
                points.append(u.cp())
                el = u.tov(v).mag()
                if el > maxlength:
                    n = 1
                    while el / n > maxlength:
                        n += 1
                    points.extend(u.line(v, n - 1))
            loops.append(points)
        return self.__class__(loops)

    def offset(self, r, closed=True, r0=10000):
        """Compute offset of self using distance r"""
        # TODO: handle holes
        # TODO: why support `closed` parameter?
        Z = vec3.Z()
        self._toxy()
        offsets = []
        for ps in self.loops:
            offset = []
            for x, y, z in slide(ps, 3):
                a = (y - x).axy(z - y)
                if isnear(a, 0, epsilon=0.01):
                    p = y + ((z - x).crs(Z).nrm() * -r)
                    offset.append(p)
                elif isnear(a, np.pi, epsilon=0.01):
                    p = y + ((z - y).crs(Z).nrm() * r)
                    offset.append(p)
                    p = y + ((y - z).crs(Z).nrm() * r)
                    offset.append(p)
                else:
                    xy = (y - x).nrm()
                    yz = (z - y).nrm()
                    u = xy.crs(Z)
                    v = yz.crs(Z)
                    a = x + (u * -r) + (xy * -r0)
                    b = y + (u * -r) + (xy *  r0)
                    c = y + (v * -r) + (yz * -r0)
                    d = z + (v * -r) + (yz *  r0)
                    p = sintsxyp(a, b, c, d)
                    offset.append(p)
            offset.insert(0, offset.pop(-1))

            if not closed:
                offset[ 0] = loop[ 0] + Z.crs(loop[ 1] - loop[ 0]).nrm() * r
                offset[-1] = loop[-1] + Z.crs(loop[-1] - loop[-2]).nrm() * r

            offsets.append(offset) #??

        self._fromxy(*offsets)
        return self.__class__(offsets)

    def union(self, other):
        """Find the union of two loops"""
        # TODO: handle holes
        assert len(self.holes) == 0
        self._toxy()
        other._toxy()
        pg = planargraph(loops=(self.loops + other.loops))
        loops = [[pg.vertices[i].cp() for i in l] for l in pg.loops()]
        areas = [loop_area(l) for l in loops]
        # each island yields exactly one loop of area < 0
        parts = [l for l, a in zip(loops, areas) if (a < 0)]
        # remove any part completely contained within another
        bad = [False] * len(parts)
        for i, u in enumerate(parts):
            for j, v in enumerate(parts):
                if not i == j:
                    if all(p.inbxy(u, False) for p in v):
                        bad[j] = True
        parts = [p for p, b in zip(parts, bad) if not b]
        self._fromxy(*parts)
        other._fromxy()
        return self.__class__(parts)#.rot(q)

    def intersect(self, other):
        """Find the intersection of two loops"""
        # TODO: handle holes
        assert len(self.holes) == 0
        self._toxy()
        other._toxy()
        pg = planargraph(loops=(self.loops + other.loops))
        _, iloops = pg.polygon()
        parts = []
        for loop in iloops:
            for p in loop:
                if self.contains(p) and other.contains(p):
                    continue
                else:
                    break
            else:
                parts.append(loop)
        self._fromxy(*parts)
        other._fromxy()
        return self.__class__(parts)#.rot(q)

    def difference(self, other, inplace=True):
        """Find the difference of two loops"""
        # TODO: handle holes
        assert len(self.holes) == 0
        self._toxy()
        other._toxy()
        pg = planargraph(loops=(self.loops + other.loops))
        _, iloops = pg.polygon()
        parts = []
        for loop in iloops:
            for p in loop:
                if not other.contains(p):
                    parts.append(loop)
                    break
        self._fromxy(*parts)
        other._fromxy()
        if inplace:
            self.loops = []
            for part in parts:
                self.al(part)
            return self#.rot(q)
        else:
            return self.__class__(parts)#.rot(q)
