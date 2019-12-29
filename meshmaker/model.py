from collections import defaultdict
from .base import Base
from .mesh import Trimesh


class Model(Base):

    def __init__(self, meshes=None, **kws):
        super().__init__(**kws)
        self.meshes = defaultdict(list) if meshes is None else meshes

    def add(self, material, mesh=None):
        mesh = Trimesh() if mesh is None else mesh
        self.meshes[material].append(mesh)
        return mesh
