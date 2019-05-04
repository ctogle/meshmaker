import numpy as np


class quat:

    def __repr__(self):
        return f'vec3({self.w:.4f}, {self.x:.4f}, {self.y:.4f}, {self.z:.4f})'

    def __init__(self, w, x, y, z):
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def av(cls, a, v):
        v = v * (np.sin(a / 2.0) / v.mag())
        return cls(np.cos(a / 2.0), v.x, v.y, v.z)

    def fp(self):
        return quat(-self.w, self.x, self.y, self.z)



