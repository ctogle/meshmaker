from .vec3 import vec3
from .quat import quat
from .geometry import loop_normal


class BSP:
    """Reimplementation of: https://github.com/evanw/csg.js/blob/master/csg.js
    Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
    by picking a polygon to split along. That polygon (and all other coplanar
    polygons) are added directly to that node and the other polygons are added to
    the front and/or back subtrees. This is not a leafy BSP tree since there is
    no distinction between internal and leaf nodes."""

    @property
    def W(self):
        return self.O.dot(self.N)

    def __init__(self, loops):
        self.O = None
        self.N = None
        self.L = None
        self.R = None
        self.loops = []
        self.build(loops)

    def cp(self):
        """Create an independent copy of this BSP instance."""
        bsp = BSP([])
        bsp.O = self.O.cp()
        bsp.N = self.N.cp()
        bsp.L = self.L.cp() if self.L is not None else None
        bsp.R = self.R.cp() if self.R is not None else None
        bsp.loops = [[p.cp() for p in l] for l in self.loops]
        return bsp

    @classmethod
    def from_mesh(self, mesh):
        """Create a BSP from a brep mesh"""
        loops = [[mesh.vertices[v].cp() for v in face] for f, face in mesh]
        return BSP(loops)

    @classmethod
    def from_plane(self, O, N, r=1000):
        """Create a BSP from a division plane"""
        loop = vec3.O().ring(r, 4)
        O.trnps(quat.toxy(N).fp().rot(loop))
        bsp = BSP([loop])
        bsp.O = O
        bsp.N = N
        return bsp

    def split_loop(self, loop, coL, coR, L, R, eps=0.000001):
        """Split `polygon` by this plane if needed, then put the polygon or polygon
        fragments in the appropriate lists. Coplanar polygons go into either
        `coL` or `coR` depending on their orientation with respect to this plane.
        Polygons in front or in back of this plane go into either `front` or `back`.
        """
        COPLANAR, FRONT, BACK, SPANNING = 0, 1, 2, 3
        ltype = 0
        types = []
        for p in loop:
            t = self.N.dot(p) - self.W
            ptype = (BACK if t < -eps else (FRONT if t > eps else COPLANAR))
            ltype |= ptype
            types.append(ptype)
        if ltype == COPLANAR:
            if self.N.dot(loop_normal(loop)) > 0:
                coL.append(loop)
            else:
                coR.append(loop)
        elif ltype == FRONT:
            L.append(loop)
        elif ltype == BACK:
            R.append(loop)
        elif ltype == SPANNING:
            l, r = [], []
            n = len(loop)
            for i in range(n):
                j = (i + 1) % n
                ti, tj = types[i], types[j]
                vi, vj =  loop[i],  loop[j]
                if ti != BACK:
                    l.append(vi)
                if ti != FRONT:
                    r.append(vi.cp() if ti != BACK else vi)
                if (ti | tj) == SPANNING:
                    t = (self.W - self.N.dot(vi)) / self.N.dot(vj - vi)
                    v = vi.lerp(vj, t)
                    l.append(v)
                    r.append(v.cp())
            if len(l) >= 3:
                L.append(l)
            if len(r) >= 3:
                R.append(r)

    def invert(self):
        """Convert solid space to empty space and empty space to solid space."""
        for loop in self.loops:
            loop.reverse()
        self.N = self.N.fp()
        if self.L:
            self.L.invert()
        if self.R:
            self.R.invert()
        self.L, self.R = self.R, self.L

    def clip_loops(self, loops):
        """Recursively remove all polygons in `polygons` that are inside this BSP tree."""
        if self.N is None:
            return loops[:]
        L, R = [], []
        for loop in loops:
            self.split_loop(loop, L, R, L, R)
        if self.L:
            L = self.L.clip_loops(L)
        if self.R:
            R = self.R.clip_loops(R)
        else:
            R = []
        return L + R

    def clip_to(self, bsp):
        """Remove all polygons in this BSP tree that are inside the other BSP tree `bsp`."""
        self.loops = bsp.clip_loops(self.loops)
        if self.L:
            self.L.clip_to(bsp)
        if self.R:
            self.R.clip_to(bsp)

    def all_loops(self):
        """Return a list of all polygons in this BSP tree."""
        loops = self.loops[:]
        if self.L:
            loops.extend(self.L.all_loops())
        if self.R:
            loops.extend(self.R.all_loops())
        return loops

    def build(self, loops):
        """Build a BSP tree out of `polygons`. When called on an existing tree, the
        new polygons are filtered down to the bottom of the tree and become new
        nodes there. Each set of polygons is partitioned using the first polygon
        (no heuristic is used to pick a good split)."""
        if not loops:
            return
        if self.N is None:
            self.O = loops[0][0].cp()
            self.N = loop_normal(loops[0])
        L, R = [], []
        for loop in loops:
            self.split_loop(loop, self.loops, self.loops, L, R)
        if L:
            if self.L is None:
                self.L = BSP([])
            self.L.build(L)
        if R:
            if self.R is None:
                self.R = BSP([])
            self.R.build(R)

    def union(self, other):
        """Add other to self"""
        self.clip_to(other)
        other.clip_to(self)
        other.invert()
        other.clip_to(self)
        other.invert()
        self.build(other.all_loops())
        return self

    def difference(self, other):
        """Subtract other from self"""
        self.invert()
        self.clip_to(other)
        other.clip_to(self)
        other.invert()
        other.clip_to(self)
        other.invert()
        self.build(other.all_loops())
        self.invert()
        return self

    def intersect(self, other):
        """Subtract the difference of self and other from self"""
        self.invert()
        other.clip_to(self)
        other.invert()
        self.clip_to(other)
        other.clip_to(self)
        self.build(other.all_loops())
        self.invert()
        return self

    def split(self, O, N):
        """Split a BSP tree into two trees using a plane"""
        L, R = self.cp(), self.cp()
        L = L.difference(BSP.from_plane(O, N))
        R = R.difference(BSP.from_plane(O.cp(), N.fp()))
        return L, R
