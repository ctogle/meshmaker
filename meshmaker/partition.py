from .base import Base
from .vec3 import vec3
from .mesh import Mesh
from .loop import Loops
from .geometry import isnear
from collections import defaultdict


class Partition(Base):

    def __init__(self, *meshes, **kws):
        super().__init__(**kws)
        self.meshcount = 0
        self.meshes = []
        for mesh in meshes:
            self.av(mesh)

    def __iter__(self):
        """Yield the index, mesh pairs in the partition"""
        for i, mesh in enumerate(self.meshes):
            if mesh is not None:
                yield (i, mesh)

    def av(self, mesh, **kws):
        """Add new volume/vertex"""
        m = len(self.meshes)
        self.meshes.append(mesh)
        self.meshcount += 1
        return m

    def rv(self, v):
        """Remove volume/vertex v"""
        mesh = self.meshes[v]
        self.meshes[v] = None
        self.meshcount -= 1
        return mesh

    def sv(self, O, N, *vs):
        """Split vertices *vs using plane defined by O and N"""
        nvs = []
        for v in vs:
            v = self.rv(v)
            x, y = v.split(O, N)
            nvs.append(self.av(x))
            nvs.append(self.av(y))
        return nvs

    def graph(self, Interface=None):
        """Map the adjacency of volumes within the partition
        identifying the intersections of their bounding faces"""
        adjacent = defaultdict(lambda : defaultdict(lambda : {}))
        normals = {i: mesh.face_normals() for i, mesh in self}
        shell = getattr(self, '_shell', None)
        for i, u in self:
            for k, uF in u:
                uN = normals[i][k]
                uL = [u.vertices[x].cp() for x in uF]
                uL.reverse()
                for j, v in self:
                    if i < j:
                        continue
                    for l, vF in v:
                        vN = normals[j][l]
                        vL = [v.vertices[x].cp() for x in vF]
                        if uN.isnear(vN.fp()):
                            if isnear(uL[0].dot(uN), vL[0].dot(uN)):
                                overlap = Loops([uL]).intersect(Loops([vL]))
                                if overlap.loops:
                                    interior = not (i == shell or j == shell)
                                    support = dict(extent=overlap, interior=interior)
                                    if Interface is not None:
                                        support = Interface(**support)
                                    adjacent[i][j][k] = support
                                    adjacent[j][i][l] = support
        return adjacent

    def shell(self):
        """Add node for the exterior of the partition (do last if at all)"""
        shell = Mesh.Union(*filter(None, self.meshes)).fp()
        self._shell = self.av(shell)
        return self._shell

    @property
    def shelled(self):
        return getattr(self, '_shell', None)
