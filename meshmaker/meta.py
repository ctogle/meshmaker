from .vec3 import vec3
from .seam import Seam
from .tform import TForm
from .model import Model
from .mesh import Mesh
from .pmesh import MetaMesh
from .geometry import near, isnear, slide, loop_normal
import numpy as np


class Stairs(MetaMesh):

    def steps(self, lb, rb):
        h = lb[1].d(lb[0])
        z = lb[1].z - lb[0].z
        t = np.arccos(z / h)
        steps = Mesh()
        steps.uvs = {}
        for (aa, bb), (cc, dd) in slide(list(zip(lb, rb)), 2, 1):
            dX = vec3.Z().crs(bb - aa).nrm() * (h * np.sin(t))
            dZ = vec3.Z(lb[1].z - lb[0].z)
            assert dd.z > aa.z, 'Stair splines must point upward!'
            front = steps.af([aa, bb, bb + vec3.Z(dd.z - bb.z), aa + vec3.Z(cc.z - aa.z)])
            top   = steps.af([aa + vec3.Z(cc.z - aa.z), bb + vec3.Z(dd.z - bb.z), dd, cc])
            lside = steps.af([aa, aa + vec3.Z(cc.z - aa.z), cc])
            rside = steps.af([dd, bb + vec3.Z(dd.z - bb.z), bb])
            seams = steps.perimeter((front, top))
            steps.unwrap_uvs(top, O=(aa + vec3.Z(cc.z - aa.z), vec3.O()), S=vec3.U(2),
                             seams=seams, uvs=steps.uvs)
        return TForm.from_meshes(steps)

    def ramp(self, lb, rb, dz=-0.1):
        dz = vec3.Z(dz)
        ramp = Mesh()
        ramp.af([lb[ 0] + dz, rb[ 0] + dz, rb[ 0], lb[ 0]])
        ramp.af([rb[-1] + dz, lb[-1] + dz, lb[-1], rb[-1]])
        for ((aa, bb), (cc, dd)) in slide(list(zip(lb, rb)), 2, 1):
            ramp.af([cc, aa, bb, dd])
            ramp.af([aa + dz, cc + dz, dd + dz, bb + dz])
            ramp.af([cc + dz, aa + dz, aa, cc])
            ramp.af([dd + dz, dd, bb, bb + dz])
        ramp.uvs = ramp.vertex_uvs(
            O=(lb[0], vec3.O()), S=vec3.U() * 2,
            seams=ramp.angle_seams())
        return TForm.from_meshes(ramp, texture='generic_13')

    def splines(self, a, b, c, d, stepheight=0.1, lmargin=0.001, rmargin=0.001):
        r = a.d(d)
        z = d.z - a.z
        n = int(abs(z) / stepheight)
        N = loop_normal([a, b, c, d])

        ldu = N.crs(b - a).nrm()
        ldv = N.crs(d - c).nrm()
        rdu = ldu.cp()
        rdv = ldv.cp()

        dw = (d - a).crs(vec3.Z()).nrm().scl(lmargin)
        lb = Seam((a + dw), (d + dw), n, ldu, ldv).loop()

        dw = (d - a).crs(vec3.Z()).nrm().scl(rmargin)
        rb = Seam((b - dw), (c - dw), n, rdu, rdv).loop()

        if z < 0:
            lb.reverse()
            rb.reverse()

        return lb, rb

    height = 0.04

    def scene(self):
        vs, fs = self.control.vertices, self.control.faces
        scene = TForm()
        for i, ((a, b), (c, d)) in enumerate(self.control._stairs):
            a, b, c, d = vs[a], vs[b], vs[c], vs[d]
            lb, rb = self.splines(a, b, d, c, stepheight=self.height)
            scene.add(self.ramp(lb, rb, dz=-self.height))
            if not (isnear(a.z, c.z) and isnear(b.z, d.z)):
                scene.add(self.steps(lb, rb))
        return scene


class Railings(MetaMesh):

    def railing(self, a, b, dz=0.4, dw=0.02):
        dz = vec3.Z(dz)
        railing = Mesh.prism([a, b, b + dz, a + dz], dw, 2)
        railing.uvs = railing.vertex_uvs(
            O=(a, vec3.O()), S=vec3.U() * 2,
            seams=railing.angle_seams())
        return TForm.from_meshes(railing, texture='generic_9')

    def scene(self):
        vs, fs = self.control.vertices, self.control.faces
        scene = TForm()
        for i, (a, b) in enumerate(self.control._railings):
            scene.add(self.railing(vs[a], vs[b]))
        return scene


'''
class Railing(MetaMesh):

    """Provide a control mesh with a `_railings` attribute"""

    def simple(self, u, v):
        x, y = u + vec3.Z(0.4), v + vec3.Z(0.4)
        mesh = Mesh.prism([u, v, y, x], 0.05, 0.0)
        mesh.uvs = mesh.vertex_uvs(
            O=(u, vec3.O()), S=vec3.U() * 2,
            seams=mesh.angle_seams())
        return TForm.from_meshes(mesh)

    def scene(self):
        chunks = []
        for u, v in self.control._railings:
            u = u if isinstance(u, vec3) else self.control.vertices[u]
            v = v if isinstance(v, vec3) else self.control.vertices[v]
            chunks.append(self.simple(u, v))
        return TForm(children=chunks)


class Stairs(MetaMesh):

    # TODO: support edge to edge staircases

    def switchback(self, u, v, z, m):
        u, v, w = u, u.lerp(v, 0.5), v

        N = vec3.Z().crs(v - u).nrm()

        dz = (u.z - z)
        dZ = vec3.Z(dz)
        dH = vec3.Z(0.1)
        dy = (u.d(w) / 2)

        o = N * (m * abs(dz / 2)) - dZ * 0.5
        p, q, r = u + o, v + o, w + o
        a, b = u.lerp(v, 0.5), p.lerp(q, 0.5)
        ll = [a, b, b - dH, a - dH]
        a, b = v.lerp(w, 0.5) - dZ, q.lerp(r, 0.5)
        rl = [a, b, b - dH, a - dH]
        ml = [q, q + N * dy, q + N * dy - dH, q - dH]
        meshes = (Mesh.prism(ll, u.d(v)),
                  Mesh.prism(rl, v.d(w)),
                  Mesh.prism(ml, u.d(w)))

        for mesh in meshes:
            mesh.uvs = mesh.vertex_uvs(
                O=(u, vec3.O()), S=vec3.U() * 2,
                seams=mesh.angle_seams())

        self.control._railings.append((p, u))
        self.control._railings.append((v, q))
        self.control._railings.append((p + N * dy, p))
        self.control._railings.append((r + N * dy, p + N * dy))
        self.control._railings.append((r, r + N * dy))
        self.control._railings.append((q, v - dZ))
        self.control._railings.append((w - dZ, r))
        self.control._railings.append((v - dZ, v - dZ - N * 0.1))
        self.control._railings.append((w - dZ - N * 0.1, w - dZ))
        return TForm.from_meshes(*meshes)

    def simple(self, u, v, p, q):
        a, b = u.lerp(v, 0.5), p.lerp(q, 0.5)
        l = [a, b, b - vec3.Z(0.1), a - vec3.Z(0.1)]
        mesh = Mesh.prism(l, u.d(v))
        mesh.uvs = mesh.vertex_uvs(
            O=(u, vec3.O()), S=vec3.U() * 2,
            seams=mesh.angle_seams())
        N = vec3.Z().crs(u - v).nrm()
        self.control._railings.append((u, p))
        self.control._railings.append((p, p + N * 0.1))
        self.control._railings.append((q, v))
        self.control._railings.append((q + N * 0.1, q))
        steps = self.steps(p, q, v, u)
        return TForm.from_meshes(mesh).add(steps).parent

    def steps(self, a, b, c, d, stepheight=0.1):
        r = a.d(d)
        z = d.z - a.z
        n = int(abs(z) / stepheight)
        lb = a.line(d, n, True)
        rb = b.line(c, n, True)
        if z < 0:
            lb.reverse()
            rb.reverse()
        h = lb[1].d(lb[0])
        z = lb[1].z - lb[0].z
        t = np.arccos(z / h)
        dX = vec3.Z().crs(b - a).nrm() * (h * np.sin(t))
        steps = Mesh()
        steps.uvs = {}
        for (aa, bb), (cc, dd) in slide(list(zip(lb, rb)), 2, 1):
            if d.z > a.z:
                O = cc - dX
                front = steps.af([aa, bb, dd - dX, cc - dX])
                top   = steps.af([cc - dX, dd - dX, dd, cc])
            else:
                O = dd + dX
                front = steps.af([bb, aa, cc + dX, dd + dX])
                top   = steps.af([dd + dX, cc + dX, cc, dd])
            seams = steps.perimeter((front, top))
            steps.unwrap_uvs(top, O=(O, vec3.O()), S=vec3.U(2),
                             seams=seams, uvs=steps.uvs)
        return TForm.from_meshes(steps)

    def place(self, u, v, w=1, align=1.0):
        if (u, v) in self.control._railings:
            self.control._railings.remove((u, v))
        u = self.control.vertices[u] if isinstance(u, int) else u
        v = self.control.vertices[v] if isinstance(v, int) else v
        d = u.d(v)
        dx = w / 2 / d
        x = near(near(align, 1, dx), 0, dx)
        if d > w:
            if x == 0:
                u_ = v - (v - u).nrm() * w
                self.control._railings.append((u, u_))
                #self.control._railings.append((u, u_.lerp(v, 0.5)))
                u = u_
            elif x == 1:
                v_ = u + (v - u).nrm() * w
                self.control._railings.append((v_, v))
                v = v_
            else:
                u_ = u.lerp(v, x - dx)
                v_ = u.lerp(v, x + dx)
                self.control._railings.append((u, u_))
                self.control._railings.append((v_, v))
                u = u_
                v = v_
        return u, v

    def scene(self):
        vs, fs = self.control.vertices, self.control.faces

        chunks = []

        w = 1.5
        x = 0.6
        m = 2.0
        #z = 0.0
        #z = vs[fs[self.patch[0]][0]].z
        sb = False

        wires = []

        for i, (tip, tail) in enumerate(self.control._stairs):
            print(tip, tail)

            #z = vs[fs[tip][0]].z
            z = vs[fs[tip[0]][0]].z
            u, v = tail

            u, v = self.place(u, v, w, x)

            wires.append([u, v])

            if sb:
                chunks.append(self.switchback(v, u, z, m))

            else:

                N = vec3.Z().crs(u - v).nrm()
                p = u + N * (m * abs(u.z - z)) - vec3.Z(u.z - z)
                q = v + N * (m * abs(v.z - z)) - vec3.Z(v.z - z)

                chunks.append(self.simple(u, v, p, q))
                #chunks.append(self.simple(v, u, z, m))

        guide = Mesh()
        guide.wires = lambda tf: [tf.transform(wire) for wire in wires]
        guide = Model(meshes={'generic_0': [guide]})

        return TForm(children=chunks, models=[guide])
'''
