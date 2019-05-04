from .vec3 import vec3
from .quat import quat
import numpy as np


class tform:

    def __repr__(self):
        p = self.position.__repr__()
        q = self.rotation.__repr__()
        s = self.scale.__repr__()
        return 'tform:\n {}\n {}\n {}'.format(p, q, s)

    def __init__(self, position=None, rotation=None, scale=None, parent=None):
        self.position = position if position else vec3(0, 0, 0)
        self.rotation = rotation if rotation else quat.av(0, vec3(0, 0, 1))
        self.scale = scale if scale else vec3(1, 1, 1)
        self.parent = parent

    def _transform(self, world):
        local = (world * self.scale).rot(self.rotation) + self.position
        if self.parent:
            return self.parent._transform(local)
        else:
            return local

    def _inverse(self, local):
        world = (local + self.position.fp()).rot(self.rotation.fp()) * self.scale.inv()
        if self.parent:
            return self.parent._inverse(world)
        else:
            return world

    def transform(self, world):
        if isinstance(world, list):
            return list(self.transform(p) for p in world)
        elif isinstance(world, tuple):
            return tuple(self.transform(p) for p in world)
        else:
            return self._transform(world)

    def inverse(self, local):
        if isinstance(local, list):
            return list(self.inverse(p) for p in local)
        elif isinstance(local, tuple):
            return tuple(self.inverse(p) for p in local)
        else:
            return self._inverse(local)


if __name__ == '__main__':
    tf = tform(
        vec3(1, 2, 0),
        quat.av(np.pi / 6, vec3(0, 0, 1)),
        vec3(2, 2, 2))
    print(tf.transform(vec3(0, 0, 0)))
    print(tf.transform(vec3(1, 1, 1)))

    root = tform(vec3(2, 1, -1), quat.av(np.pi / 6, vec3(0, 0, 1)), vec3(1, 2, 2))
    node = tform(vec3(0, 0, 5), quat.av(np.pi / 3, vec3(0, 0, 1)), vec3(1, 1, 1), root)
    print(root.transform(vec3(0, 0, 0)))
    print(node.transform(vec3(0, 0, 0)))
    print(root.transform(vec3(1, 1, 1)))
    print(node.transform(vec3(1, 1, 1)))
