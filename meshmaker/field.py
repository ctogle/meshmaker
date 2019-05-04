
from dilap.geometry import vec3
from dilap.core.plotting import *
import numpy

class basis(vec3):

    epsilon = 0.00000001

    @classmethod
    def from_r_theta(cls, r, th):
        a = r * numpy.cos(2 * th)
        b = r * numpy.sin(2 * th)
        return cls(a, b, 0)

    @classmethod
    def from_ab(cls, a, b):
        return cls(a, b, 0)

    def evalues(self):
        d = self.mag()
        return -d, d

    def evectors(self):
        if abs(self.y) < basis.epsilon:
            if abs(self.x) < basis.epsilon:
                #major, minor = vec3(0, 0, 0), vec3(0, 0, 0)
                major, minor = None, None
            else:
                major, minor = vec3(1, 0, 0), vec3(0, 1, 0)
        else:
            e1, e2 = self.evalues()
            major = vec3(self.y, e1 - self.x, 0)
            minor = vec3(self.y, e2 - self.x, 0)
        return major, minor

    def plot(self, p, ax, which=None):
        l = 1.5
        major, minor = self.evectors()
        if (which == 'major' or not which) and major:
            major.nrm()
            dp = major.cp().uscl(l)
            plot_point_xy(p + dp, ax, col='b')
            plot_edges_xy((p, p + dp), ax, col='b', lw=2)
        if (which == 'minor' or not which) and minor:
            minor.nrm()
            dp = minor.cp().uscl(l)
            plot_point_xy(p + dp, ax, col='g')
            plot_edges_xy((p, p + dp), ax, col='g', lw=2)
        plot_point_xy(p, ax, col='r')
		return ax


from .vec3 import vec3
from .basis import basis
import numpy

def field_element(p_0, decay, weight):
    def wrap(f):
        def wrapped(p):
            return f(p).uscl(weight * numpy.exp(-decay * p_0.dxy(p) ** 2))
        return wrapped
    return wrap

def radial_element(p_0, decay, weight):
    @field_element(p_0, decay, weight)
    def getbasis(p):
        dpx = (p.x - p_0.x)
        dpy = (p.y - p_0.y)
        a = dpy * dpy - dpx * dpx
        b = - 2 * dpx * dpy
        return basis.from_ab(a, b).nrm()
    return getbasis

def grid_element(p_0, decay, weight, major):
    @field_element(p_0, decay, weight)
    def getbasis(p):
        d = major.cp().nrm()
        th = numpy.arctan2(d.y, d.x) + numpy.pi / 2.0
        return basis.from_r_theta(1, th).nrm()
    return getbasis

def topographical_element(weight, image, tform):
    xgrad, ygrad = numpy.gradient(image)
    ry, rx = image.shape
    @field_element(vec3(0, 0, 0), 0.0, weight)
    def getbasis(wp):
        ix, iy, iz = tform.wtoi(*wp)
        i, j = min(rx - 1, int(rx * ix)), min(ry - 1, int(ry * iy))
        gx, gy = xgrad[j, i], ygrad[j, i]
        r = numpy.sqrt(gx ** 2 + gy ** 2)
        th = numpy.arctan2(gy, gx) + numpy.pi / 2.0
        return basis.from_r_theta(r, th).nrm()
    return getbasis

def parse_element(element, decay, weight):
    if isinstance(element, tuple):
        yield topographical_element(weight, *element)
    elif isinstance(element, list):
        for j in range(1, len(element)):
            u, v = element[j - 1], element[j]
            yield grid_element(u.mid(v), decay, weight, u.tov(v).nrm())
    else:
        yield radial_element(element, decay, weight)

def parse_elements(elements):
    field = []
    for e, d, w in elements:
        field.extend(list(parse_element(e, d, w)))
    return field

class tensor_field:

    def __init__(self, elements):
        self.elements = parse_elements(elements)
        self.layers = len(self.elements)

    def basis(self, p):
        t = self.elements[0](p)
        if self.layers > 1:
            for e in self.elements[1:]:
                t.trn(e(p))
        return t.nrm()

    def trace(self, seed, max_length, max_steps, mode='major'):

        def orient():
            t = self.basis(p)
            major, minor = t.evectors()
            if mode == 'major':
                return major
            elif mode == 'minor':
                return minor
            else:
                raise ValueError('unknown mode: "%s"' % mode)

        def forward(last_dp):
            direction = orient()
            if direction:
                if last_dp.dot(direction) < 0:
                    direction.flp()
                dp = direction.nrm().uscl(max_length)
            else:
                dp = vec3(0, 0, 0)
            p.trn(dp)
            return dp

        p, dp = seed
        steps = 0
        while steps < max_steps:
            yield p.cp(), dp.cp()
            dp = forward(dp)
            steps += 1

    def plot(self, ax, ps, which=None):
        for p in ps:
            t = self.basis(p)
            t.plot(p, ax, which=which)
        return ax


