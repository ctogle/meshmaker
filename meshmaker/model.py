from collections import defaultdict
from .base import Base
from .mesh import Mesh


class Model(Base):

    @classmethod
    def cube_model(cls, material, r=1):
        mesh = Mesh.cube_mesh(r=r)
        model = cls()
        model.add(material, mesh)
        return model

    def __init__(self, meshes=None, **kws):
        super().__init__(**kws)
        self.meshes = defaultdict(list) if meshes is None else meshes

    def add(self, material, mesh=None):
        mesh = Mesh() if mesh is None else mesh
        self.meshes[material].append(mesh)
        return mesh
