from meshmaker.geometry import slide


class Seam:
    """Parameterized curve (optionally a catmull-rom spline)"""

    @classmethod
    def l2s(cls, loop, m=0):
        s = [cls(u, v) for u, v in slide(loop, 2, m)]
        for i in range(1, len(s)):
            s[i - 1].next = s[i]
        return s[0]

    def slide(self, n=2, m=1):
        # this doesnt work because __iter__ cant receive n
        yield from slide(list(self), n=n, m=m)

    def __init__(self, u, v, n=0, du=None, dv=None):
        self.u = u
        self.v = v
        self.n = n
        self.du = (v - u).nrm() if du is None else du.cp()
        self.dv = (u - v).nrm() if dv is None else dv.cp()

    def __iter__(self):
        yield from self._iter()

    def _iter(self, n=None):
        m = self.n if n is None else n
        if m == 0:
            yield self.u
            #yield self.v
        else:
            yield self.u
            yield from self.u.spline(self.v, self.du, self.dv, m, alpha=1)
        if not (hasattr(self, 'next') and self.next.u.isnear(self.v)):
            yield self.v
        if hasattr(self, 'next'):
            yield from self.next._iter(n=n)

    def _control(self):
        unique = []
        for p in self._iter(n=0):
            if not p in unique:
                unique.append(p)
        return unique

    def trn(self, t):
        t.trnps(self._control())
        return self

    def scl(self, s):
        s.sclps(self._control())
        return self

    def rot(self, q):
        q.rot(self._control())
        return self

    def cp(self, root=True):
        cp = Seam(self.u.cp(), self.v.cp(), self.n,
                  self.du.cp(), self.dv.cp())
        if hasattr(self, 'next'):
            cp.next = self.next.cp(False)
            cp.next.u = cp.v
        if root and self._closed():
            cp._last().v = cp.u
        return cp

    def _last(self):
        return self.next._last() if hasattr(self, 'next') else self

    def _closed(self):
        return self._last().v.isnear(self.u)

    def loop(self, n=None):
        points = list(self._iter(n=n))
        if self._closed():
            points.pop(-1)
        return points
