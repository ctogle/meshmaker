from .base import Base
from .mat44 import mat44


class TForm(Base):

    def world(self):
        # NOTE: TForms cant be reused because they can have exactly one parent
        if self.parent is None:
            base = mat44.I()
        else:
            base = self.parent.world()
        return base * self.transformation

    def translate(self, t):
        self.transformation.m14 += t.x
        self.transformation.m24 += t.y
        self.transformation.m34 += t.z

    def transform(self, vs):
        w = self.world()
        if isinstance(vs, (list, tuple)):
            return [w * v for v in vs]
        else:
            return w * vs

    def __repr__(self):
        return self.transformation.__repr__()

    def __init__(self, t=None, q=None, s=None,
                 parent=None, children=None, **kws):
        super().__init__(**kws)
        self.transformation = mat44.tform(t, q, s)

        self.parent = parent
        if parent:
            parent.add(self)
        self.children = []
        if children is not None:
            for child in children:
                self.add(child)

    def cp(self):
        cp = TForm()
        for k, v in self.__dict__.items():
            if k == 'children':
                for c in v:
                    cp.add(c.cp())
            else:
                setattr(cp, k, v)
        return cp

    def add(self, child):
        child.parent = self
        if not child in self.children:
            self.children.append(child)
        return child
