from .vec3 import vec3
from .quat import quat
from .geometry import slide
from .planargraph import planargraph
from .delaunay import triangulation
from contextlib import contextmanager


class loop:
    # TODO: lazy/refreshable interface??

    @property
    def N(self):
        """Compute the normal vector of self"""
        N = vec3.O()
        for u, v, w in slide(self.ps, 3):
            uv, vw = (v - u), (w - v)
            N.trn(uv.crs(vw))
        assert not N.isO(), 'zero vector in loop.N ...'
        return N.nrm()

    @contextmanager
    def orientation(self, N):
        """For operations in a particular plane..."""
        q = quat.toxy(self.N)
        self.rot(q)
        yield q.fp()
        self.rot(q.fp())

    @property
    def area(self):
        """Compute the area contained by self"""
        with self.orientation(vec3.Z()):
            area = 0.0
            for u, v in slide(self.ps, 2):
                area -= (u.x + v.x) * (u.y - v.y) / 2.0
        return area

    #def containstri(self, a, b, c):
    #    # query
    #    raise

    def triangulation(self, e=0.0001, h=None, r=10000):
        """Compute a delaunay triangulation of self"""
        # query
        with self.orientation(vec3.Z()) as q:
            py = (self.ps, [h.ps for h in self.holes])
            t = triangulation(py, e, h, r)
            q.rot(t.points)
        return t

    def __len__(self):
        return self.len

    def __getitem__(self, i):
        return self.ps[i]

    def __iter__(self):
        return self.ps.__iter__()

    def __init__(self, points, holes=None):
        self.ps = points
        self.len = len(self.ps)
        self.holes = [] if holes is None else holes

    def rot(self, q):
        q.rot(self.ps)
        for hole in self.holes:
            hole.rot(q)
        return self

    def add(self, hole):
        self.holes.append(hole)
        return self

    def cp(self):
        points = [p.cp() for p in self.ps]
        holes = [h.cp() for h in self.holes]
        return self.__class__(points, holes=holes)

    def fp(self):
        self.ps.reverse()
        for h in self.holes:
            h.fp()
        return self

    def smooth(self, weight=0.5, iterations=1):
        """Smooth towards 1-ring neighbors in-place"""
        # TODO: handle holes
        for j in range(iterations):
            dps = []
            for u, v, w in slide(self.ps, 3):
                dp = (vec3.com((u, w)) - v) * weight
                dps.append(dp)
            dps.insert(0, dps.pop(-1))
            for p, dp in zip(self.ps, dps):
                p.trn(dp)
        return self

    def edgesplit(self, maxlength):
        """Compute new loop where no edge exceeds maxlength"""
        # TODO: handle holes
        points = []
        for u, v in slide(self.ps, 2):
            points.append(u.cp())
            el = u.tov(v).mag()
            if el > maxlength:
                n = 1
                while el / n > maxlength:
                    n += 1
                points.extend(u.line(v, n - 1))
        return self.__class__(points)

    def offset(self, r, closed=True, r0=10000):
        """Compute offset of self using distance r"""
        # TODO: handle holes
        # TODO: why support `closed` parameter?
        Z = vec3.Z()
        with self.orientation(Z) as q:
            offset = []
            for x, y, z in slide(self.ps, 3):
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

        return self.__class__(offset).rot(q)

    def union(self, other):
        """Find the union of two loops"""
        # TODO: handle holes
        Z = vec3.Z()
        with self.orientation(Z) as q, other.orientation(Z):
            l1 = [p.quant() for p in self.cp()]
            l2 = [p.quant() for p in other.cp()]
            segs = list(slide(l1, 2)) + list(slide(l2, 2))
            pg = planargraph(segs=segs)
            part, _ = pg.polygon()
        parts = [self.__class__(part).rot(q)]
        return parts

    def intersect(self, other):
        """Find the intersection of two loops"""
        # TODO: handle holes
        Z = vec3.Z()
        with self.orientation(Z) as q, other.orientation(Z):
            l1 = [p.quant() for p in self.cp()]
            l2 = [p.quant() for p in other.cp()]
            segs = list(slide(l1, 2)) + list(slide(l2, 2))
            pg = planargraph(segs=segs)
            _, iloops = pg.polygon()
            parts = []
            for loop in iloops:
                for p in loop:
                    if (p.inbxy(l1, True) and p.inbxy(l2, True)):
                        continue
                    else:
                        break
                else:
                    parts.append(loop)
        parts = [self.__class__(part).rot(q) for part in parts]
        return parts

    def difference(self, other):
        """Find the difference of two loops"""
        # TODO: handle holes
        Z = vec3.Z()
        with self.orientation(Z) as q, other.orientation(Z):
            l1 = [p.quant() for p in self.cp()]
            l2 = [p.quant() for p in other.cp()]
            segs = list(slide(l1, 2)) + list(slide(l2, 2))
            pg = planargraph(segs=segs)
            _, iloops = pg.polygon()
            parts = []
            for loop in iloops:
                for p in loop:
                    if not p.inbxy(l2, True):
                        parts.append(loop)
                        break
        parts = [self.__class__(part).rot(q) for part in parts]
        return parts
