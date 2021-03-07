from .base import Base
from .vec3 import vec3
from .tform import TForm
from .model import Model
from .mesh import Mesh
from .geometry import loop_normal
from collections import defaultdict
import json


def factoryfeature(f):
    """Decorator to cache function results based on kws"""
    lookup = {}
    def wrap(*ags, **kws):
        key = json.dumps(kws, sort_keys=True)
        feature = lookup.get(key)
        if feature is None:
            feature = f(*ags, **kws)
            lookup[key] = feature
        if callable(feature):
            feature.__name__ += f'_{key}'
        return feature
    return wrap


class ParamMesh(Base):
    """Proxy class for representing hierarchical
    mesh parameterizations via control meshes"""

    def __init__(self, control, **kws):
        super().__init__(features=[], **kws)
        self.control = control
        self.showcontrol = False
        self.selected = [f for f, face in control]

        from pyglet.window import key
        self.parameters = {
            key.C: self.on_toggle, key.SPACE: self.on_cycle,
            key.UP: self.on_up, key.DOWN: self.on_down,
            key.LEFT: self.on_left, key.RIGHT: self.on_right,
            key.J: self.on_j, key.K: self.on_k,
        }

    # methods for modifying control meshes
    def on_toggle(self, *ags):
        self.showcontrol = not self.showcontrol
    def on_cycle(self, *ags):
        self.selected.append(self.selected.pop(0))
    def on_left(self, *ags):
        vec3(0.5, 1, 1).sclps(self.control.vertices)
    def on_right(self, *ags):
        vec3(2, 1, 1).sclps(self.control.vertices)
    def on_up(self, *ags):
        vec3(1, 2, 1).sclps(self.control.vertices)
    def on_down(self, *ags):
        vec3(1, 0.5, 1).sclps(self.control.vertices)
    def on_j(self, *ags):
        vec3(1, 1, 0.5).sclps(self.control.vertices)
    def on_k(self, *ags):
        vec3(1, 1, 2).sclps(self.control.vertices)

    def scene(self):
        """Can be rendered as the control mesh or its parameterization.

        Parameterization is stored in list of tuples self.features, where
        entry (f, target) refers to a function `f` and some target geometry
        `target`, associated with the control mesh - f must receive the mesh
        and target and return a TForm (i.e. node in a scenegraph).
        """
        if self.showcontrol:
            s = self.selected[0]
            sface = self.control.faces[s]
            sface = [self.control.vertices[v].cp() for v in sface]
            sface = (self.control.face_normals()[s] * 0.001).trnps(sface)
            selected = Mesh()
            selected.af(sface)
            meshes = {'generic_1': [self.control],
                      'generic_8': [selected]}
            parts = [TForm(models=[Model(meshes=meshes)])]
        else:
            parts = [f(self.control, target) for f, target in self.features]
        return TForm(children=parts)

    #@factoryfeature
    @staticmethod
    def textured(texture='generic_0'):
        # TODO: raise error if texture is not a kw because of factoryfeature?
        def facefeature(control, faces):
            mesh = Mesh()
            for f in faces:
                loop = [control.vertices[v] for v in control.faces[f]]
                mesh.af(loop)
            mesh.uvs = mesh.unwrap_uvs(seams=mesh.angle_seams())
            return TForm(models=[Model(meshes={texture: [mesh]})])
        facefeature.__name__ = texture
        return facefeature


class MetaMesh(ParamMesh):
    """A ParamMesh where features are inferred from meta data
    associated with the simplices of the control mesh."""

    def metafeatures(self, **kws):
        """Aggregate feature data based on control mesh meta data."""
        features = defaultdict(list)
        methods = {}
        for f, face in self.control:
            metas = self.control.meta.get(f)
            if not isinstance(metas, tuple):
                metas = (metas, )
            for meta in metas:
                if isinstance(meta, str):
                    meta = self.textured(texture=meta)
                if callable(meta):
                    methods[meta.__name__] = meta
                    features[meta.__name__].append(f)
                elif meta is not None:
                    raise ValueError(f'whats goin on here: {meta}')
        return [(methods[m], faces) for m, faces in features.items()]

    def parameterize(self, texture='generic_0', **kws):
        """Trivial parameterization covers all faces with generic_0.
        Subclasses should add complexity."""
        generic = self.textured(texture=texture)
        for f, face in self.control:
            if self.control.meta.get(f, None) is None:
                self.control.meta[f] = generic

    #@classmethod
    #def prism(cls, loop, depth):
    #    """Loop is the projection of the extruded loop in the symmetry plane
    #    of the resulting solid, which is a prism by construction"""
    #    return cls(Mesh.prism(loop, depth))

    def __init__(self, control, **kws):
        super().__init__(control, **kws)
        self.parameterize(**kws)
        self.features = self.metafeatures(**kws)


class WireMesh(MetaMesh):
    """Convenient MetaMesh subclass for visualizing just the edges"""

    def scene(self):
        wires = []
        for i, j in self.control.e2f:
            u, v = self.control.vertices[i], self.control.vertices[j]
            wires.append((u, v))
        mesh = Mesh()
        mesh.wires = lambda tf: [tf.transform(wire) for wire in wires]
        return TForm.from_meshes(mesh)


class MetaScene(Base):
    """Hierarchy of MetaMeshes - Allows rebuilding a scene using reference
    to scenegraph (TForm root) where each node may have `metas` attribute,
    a list of MetaMesh instances"""

    def __init__(self, root, **kws):
        super().__init__(root=root, **kws)

        self.showcontrols = False

        from pyglet.window import key
        self.parameters = {
            key.C: self.on_C,
            key.I: self.on_I,
            key.J: self.on_J,
            key.K: self.on_K,
        }

    # methods which modify the entire scenegraph
    def on_C(self, key, action, modifiers):
        self.showcontrols = not self.showcontrols
    def on_I(self, key, action, modifiers):
        if modifiers.shift:
            x = vec3((0.9 if modifiers.ctrl else 1.1), 1, 1)
            self.root.scale(x)
        else:
            x = vec3.X(-1 if modifiers.ctrl else 1)
            self.root.translate(x)
    def on_J(self, key, action, modifiers):
        if modifiers.shift:
            y = vec3(1, (0.9 if modifiers.ctrl else 1.1), 1)
            self.root.scale(y)
        else:
            y = vec3.Y(-1 if modifiers.ctrl else 1)
            self.root.translate(y)
    def on_K(self, key, action, modifiers):
        if modifiers.shift:
            z = vec3(1, 1, (0.9 if modifiers.ctrl else 1.1))
            self.root.scale(z)
        else:
            z = vec3.Z(-1 if modifiers.ctrl else 1)
            self.root.translate(z)

    def build(self, cf):
        """Recursively build a scene by calling `scene` method of each
        MetaMesh instance found in the scenegraph at/below node `cf`"""
        tf = TForm.from_mat44(cf.transformation,
                              models=getattr(cf, 'models', []))
        if getattr(cf, 'metas', None) is not None:
            for meta in cf.metas:
                meta.showcontrol = self.showcontrols
                tf.add(meta.scene())
        for child in cf.children:
            tf.add(self.build(child))
        return tf

    def scene(self):
        """Recursively rebuild the scene from the root scenegraph node"""
        return self.build(self.root)
