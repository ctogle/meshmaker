from .base import Base, Laziness
from .vec3 import vec3
from .quat import quat
from .geometry import isnear, near, slide, batch, loop_normal
from .delaunay import triangulation
from .bsp import BSP
from collections import defaultdict
from functools import reduce
import numpy as np


class MeshAttribute(Laziness):
    """Lazy lookup of mesh related data"""

    def method(self, mesh):
        raise NotImplementedError

    def __init__(self, mesh, **kws):
        super().__init__(self.method(mesh), **kws)


class FaceNormals(MeshAttribute):
    """Lazy lookup of face normals"""

    @staticmethod
    def method(mesh):
        def wrapped(f):
            vs = [mesh.vertices[v] for v in mesh.faces[f]]
            u, v, w = vs[0], vs[1], vs[2]
            N = (v - u).crs(w - u).nrm()
            return N
        return wrapped


class FaceTangents(MeshAttribute):
    """Lazy lookup of face tangents"""

    @staticmethod
    def method(mesh):
        def wrapped(f):
            N = mesh.face_normals()[f]
            T = vec3.X() if isnear(abs(N.z), 1) else N.crs(vec3.Z()).nrm()
            return T
        return wrapped


class Mesh(Base):
    """2-manifold tri/quad mesh object"""

    def wires(self, tf):
        """Default to no wires"""
        return []

    def fwires(self, tf):
        """True face edges"""
        wires = []
        for f, face in self:
            ps = [self.vertices[v] for v in face]
            ps = tf.transform(ps)
            wires.extend(list(slide(ps, 2)))
        return wires

    def T2F_N3F_V3F(self, tf, smooth=False):
        """Compute a T2F_N3F_V3F representation for rendering"""
        #self.dissolve()
        normals = self.vertex_normals(smooth=smooth)
        uvs = self.vertex_uvs()
        wireframe = []
        data = []
        for f, face in self:
            ps, ns, us = [self.vertices[v] for v in face], normals[f], uvs[f][:]
            ps, ns = tf.transform(ps), tf.transform(ns)
            if len(ps) == 3:
                pass
            elif len(ps) == 4:
                ps.insert(3, ps[2]);ps.insert(3, ps[0])
                ns.insert(3, ns[2]);ns.insert(3, ns[0])
                us.insert(3, us[2]);us.insert(3, us[0])
                pass
            else:
                print(len(ps), ps)
                #raise NotImplementedError
                continue
            for p, n, u in zip(ps, ns, us):
                data.extend([u.x, u.y, n.x, n.y, n.z, p.x, p.y, p.z])

            for u, v, w in batch(ps, 3):
                wireframe.append([[u.x, u.y, u.z],
                                  [v.x, v.y, v.z],
                                  [w.x, w.y, w.z],
                                  [u.x, u.y, u.z]])
                #wireframe.extend([(u, v), (v, w), (w, u)])

        #return data, wireframe
        return data

    @classmethod
    def lattice(cls, nx=3, ny=3, nz=3):
        mesh = cls()
        X = np.linspace(-1, 1, nx) if isinstance(nx, int) else nx
        Y = np.linspace(-1, 1, ny) if isinstance(ny, int) else ny
        Z = np.linspace(-1, 1, nz) if isinstance(nz, int) else nz
        for i in range(1, len(X)):
            x0, x1 = X[i - 1], X[i]
            for j in range(1, len(Y)):
                y0, y1 = Y[j - 1], Y[j]
                for k in range(1, len(Z)):
                    z0, z1 = Z[k - 1], Z[k]
                    a, b = vec3(x0, y0, z0), vec3(x1, y0, z0)
                    c, d = vec3(x1, y1, z0), vec3(x0, y1, z0)
                    e, f = vec3(x0, y0, z1), vec3(x1, y0, z1)
                    g, h = vec3(x1, y1, z1), vec3(x0, y1, z1)
                    mesh.af([a, b, c, d], raise_on_duplicate=False)
                    mesh.af([b, a, e, f], raise_on_duplicate=False)
                    mesh.af([a, d, h, e], raise_on_duplicate=False)
                    mesh.af([e, f, g, h], raise_on_duplicate=False)
                    mesh.af([c, d, h, g], raise_on_duplicate=False)
                    mesh.af([b, c, g, f], raise_on_duplicate=False)
        mesh._X = X
        mesh._Y = Y
        mesh._Z = Z
        return mesh

    @classmethod
    def cube_mesh(cls, r=1, meta='generic_0'):
        """Generate a generate cube mesh instance"""
        mesh = cls()
        bottom = [vec3(-r,-r,-r), vec3( r,-r,-r), vec3( r, r,-r), vec3(-r, r,-r)]
        top    = [vec3(-r,-r, r), vec3( r,-r, r), vec3( r, r, r), vec3(-r, r, r)]
        mesh.bridge(bottom, top, meta=meta)
        mesh.af(top, meta=meta)
        mesh.af(bottom[::-1], meta=meta)
        return mesh

    @classmethod
    def cube_mesh2(cls, r=1, meta='generic_0'):
        mesh = cls()
        bottom = [vec3(-r,-r,-r), vec3( r,-r,-r), vec3( r, r,-r), vec3(-r, r,-r)]
        top    = [vec3(-r,-r, r), vec3( r,-r, r), vec3( r, r, r), vec3(-r, r, r)]
        mesh.fan(vec3.com(bottom), bottom[::-1], meta=meta)
        mesh.fan(vec3.com(top), top, meta=meta)
        for (u, v), (x, y) in zip(slide(bottom, 2), slide(top, 2)):
            loop = [u, v, y, x]
            mesh.fan(vec3.com(loop), loop, meta=meta)
        return mesh

    @classmethod
    def cylinder_mesh(cls, h=2, r=1, n=8, closed=False, meta='generic_0'):
        mesh = cls()
        a = vec3.Z(-h / 2)
        b = a.ring(r, n, False)
        c = a + vec3.Z(h)
        d = c.ring(r, n, False)
        mesh.bridge(b, d, meta=meta)
        if closed:
            mesh.fan(a, b[::-1], meta=meta)
            mesh.fan(c, d, meta=meta)
        return mesh

    @classmethod
    def prism(cls, loop, depth, x_align=1.0):
        """Loop is the projection of the extruded loop in the symmetry plane
        of the resulting solid, which is a prism by construction"""
        N = loop_normal(loop)
        l = (N * (x_align * depth / 2)).trnps([p.cp() for p in loop])
        r = (N * ((2 - x_align) * -depth / 2)).trnps([p.cp() for p in loop])
        C = cls()
        C.apy((l, ()))
        C.apy((r[::-1], ()))
        C.bridge(r, l)
        return C

    def cp(self):
        """Make a copy of self"""
        o = self.__class__()
        o.vertices = [p.cp() for p in self.vertices]
        o.faces = [(face[:] if face else None) for face in self.faces]
        o.nface = self.nface
        for e in self.e2f:
            o.e2f[e] = self.e2f[e]
        for v in self.v2f:
            o.v2f[v] = self.v2f[v].copy()
        for f in self.meta:
            o.meta[f] = self.meta[f]
        return o

    def fp(self):
        for f, face in self:
            face.reverse()
        return self

    def subdivide(self):
        """Create a subdivided mesh from self (currently only splits tris/quads)"""
        new = self.__class__()
        new.vertices = self.vertices[:]
        min_index = len(self.vertices)
        midpoints = {}
        for i, j in self.e2f:
            if (j, i) in midpoints:
                midpoints[(i, j)] = midpoints[(j, i)]
            else:
                midpoint = new.vertices[i].lerp(new.vertices[j], 0.5)
                midpoints[(i, j)] = new.av(midpoint, e=None, min_index=min_index)
        for f, face in self:
            if len(face) == 3:
                u, v, w = face
                p, q, r = midpoints[(u, v)], midpoints[(v, w)], midpoints[(w, u)]
                new.af([p, q, r])
                new.af([u, p, r])
                new.af([v, q, p])
                new.af([w, r, q])
            elif len(face) == 4:
                u, v, w, x = face
                p, q, r, s = (midpoints[(u, v)], midpoints[(v, w)],
                              midpoints[(w, x)], midpoints[(x, u)])
                c = vec3.com((new.vertices[u], new.vertices[v],
                              new.vertices[w], new.vertices[x]))
                c = new.av(c, e=None, min_index=min_index)
                new.af([u, p, c, s])
                new.af([p, v, q, c])
                new.af([c, q, w, r])
                new.af([s, c, r, x])
            else:
                new.af(face)
        return new

    @classmethod
    def from_bsp(cls, bsp):
        mesh = cls()
        for loop in bsp.all_loops():
            mesh.af(loop)
        return mesh

    @classmethod
    def Union(cls, *meshes):
        bsps = [BSP.from_mesh(mesh) for mesh in meshes]
        return cls.from_bsp(reduce(lambda x, y: x.union(y), bsps))

    def union(self, other):
        union = BSP.from_mesh(self).union(BSP.from_mesh(other))
        return self.__class__.from_bsp(union)

    def intersect(self, other):
        union = BSP.from_mesh(self).intersect(BSP.from_mesh(other))
        return self.__class__.from_bsp(union)

    def difference(self, other):
        union = BSP.from_mesh(self).difference(BSP.from_mesh(other))
        return self.__class__.from_bsp(union)

    def split(self, O, N):
        """Split via a plane"""
        u, v = BSP.from_mesh(self).split(O, N)
        return self.__class__.from_bsp(u), self.__class__.from_bsp(v)

    def laplacian(self):
        """Compute a laplacian matrix for self"""
        N = len(self.vertices)
        W = np.zeros((N, N))
        for i, j in self.e2f:
            W[i, j] = 1
        L = np.eye(N)
        for i, j in self.e2f:
            L[i, j] = -W[i, j] / W[i, :].sum()
        return L

    def to_xyz(self):
        """Return a vertices/faces representation"""
        xyz = np.array([[p.x, p.y, p.z] for p in self.vertices])
        faces = list(self)
        return xyz, faces

    @classmethod
    def from_xyz(cls, xyz, faces):
        return cls()._from_xyz(xyz, faces)

    def _from_xyz(self, xyz, faces):
        """Recreate Mesh in-place from vertices/faces

        Args:
            xyz (iterable): Iterable of (x, y, z) tuples
            faces (iterable): Iterable of (index, face) tuples

        """
        self.clear()
        mapping = {}
        for i, (x, y, z) in enumerate(xyz):
            mapping[i] = self.av(vec3(x, y, z))
        for f, face in faces:
            face = [mapping[v] for v in face]
            self.af(face)
        return self

    def smooth(self, inplace=False):
        """Create new Mesh with smoothing applied"""
        L = self.laplacian()
        xyz, faces = self.to_xyz()
        xyz -= np.matmul(L, xyz)
        if inplace:
            return self._from_xyz(xyz, faces)
        else:
            return self.__class__.from_xyz(xyz, faces)

    def deform(self, constraints, inplace=False):
        """Create deformation Mesh via constraints

        Args:
            constraints (dict): Dictionary mapping vertex
            index to constrained location after deformation

        Returns:
            New and deformed Mesh instance

        """
        # convert mesh to laplacian coordinates
        (xyz, faces), L = self.to_xyz(), self.laplacian()
        delta = np.matmul(L, xyz)
        # apply constraints on L/delta
        I = np.eye(len(self.vertices))
        for j, anchor in constraints.items():
            L[j, :] = I[j]
            delta[j, :] = (anchor.x, anchor.y, anchor.z)
        # recover cartesian coordinates from laplacian coordinates
        xyz, _, _, _ = np.linalg.lstsq(L, delta)
        if inplace:
            #for i, (x, y, z) in enumerate(xyz):
            #    self.vertices[i].set(x, y, z)
            #return self
            return self._from_xyz(xyz, faces)
        else:
            return self.__class__.from_xyz(xyz, faces)

    def clear(self):
        """Reinitialize data structures"""
        self.vertices = []
        self.faces = []
        self.nface = 0
        self.e2f = {}
        self.v2f = defaultdict(set)
        self.meta = {}

    def __init__(self):
        self.clear()

    def __iter__(self):
        """Yield existing faces and their indices"""
        for f, face in enumerate(self.faces):
            if face is not None:
                yield f, face

    def face_rings(self):
        """Generate mapping of faces to their neighbors (1-ring)"""
        f2f = defaultdict(list)
        for f, face in self:
            for u, v in slide(face, 2):
                adj = self.e2f.get((v, u))
                if adj is not None:
                    f2f[f].append(adj)
        return f2f

    def face_normals(self):
        """Generate face normal vector lookup"""
        if getattr(self, '_face_normals', None) is None:
            self._face_normals = FaceNormals(self)
        return self._face_normals

    def face_tangents(self, normals=None):
        """Generate face tangent vector lookup"""
        if getattr(self, '_face_tangents', None) is None:
            self._face_tangents = FaceTangents(self)
        return self._face_tangents

    def vertex_rings(self):
        """Generate mapping of vertices to their neighbors (1-ring)"""
        v2v = defaultdict(set)
        for i, j in self.e2f:
            v2v[i].add(j)
            v2v[j].add(i)
        return v2v

    def vertex_normals(self, smooth=False):
        """Generate vertex normal vector lookup with optional smoothing"""
        normals = self.normals if getattr(self, 'normals', None) else {}
        self._face_normals = None
        face_normals = self.face_normals()
        for f, face in self:
            if normals.get(f) is None:
                if smooth:
                    ns = [[face_normals[o] for o in self.v2f[v]] for v in face]
                    ns = [vec3.sum(n).nrm() for n in ns]
                else:
                    ns = [face_normals[f]] * len(face)
                normals[f] = ns
        return normals

    def vertex_uvs(self, **kws):
        """Generate UV coordinate lookup possibly unwrapping missing faces"""
        uvs = self.uvs if getattr(self, 'uvs', None) else {}
        missing = []
        for f, face in self:
            if uvs.get(f) is None:
                missing.append(f)
        if missing:
            seams = kws.pop('seams', None)
            if seams is None:
                seams = self.perimeter(missing)
                seams = set([x for y in seams for x in y])
            while missing:
                f = missing.pop(0)
                if uvs.get(f) is None:
                    self.unwrap_uvs(f, seams=seams, uvs=uvs, **kws)
        return uvs

    def unwrap_uvs(self, f=None, O=None, X=None, Y=None, S=None, seams=None, uvs=None):
        """Recursively generate UV coordinate lookup via angle based flattening
        NOTE: For sufficiently large meshes, expect recursion errors
        TODO: Reimplement this without recursion...

        Args:
            f (int): Current face index being unwrapped
            O (tuple): Origin pair (R3, UV)
            X (vec3): U-basis vector of parameterization
            Y (vec3): V-basis vector of parameterization
            S (vec3): Scale vector of parameterization
            seams (set): Set of directed edges which may not be traversed during unwrapping
            uvs (dict): Current UV mapping being computing

        Returns:
            dict mapping face indices to sequences of UV coordinates

        """
        if f is None:
            for f, face in self:
                return self.unwrap_uvs(f, O, X, Y, S, seams, uvs)
        face_normals = self.face_normals()
        N = face_normals[f]
        O = (vec3.O(), vec3.O()) if O is None else O
        if X is None and Y is None:
            X = self.face_tangents()[f]
            Y = N.crs(X)
        elif X is None:
            X = N.crs(Y) # ??? might be flipping convention
        elif Y is None:
            Y = N.crs(X) # ??? might be flipping convention
        S = vec3.U() if S is None else S
        seams = set() if seams is None else seams
        uvs = {} if uvs is None else uvs
        locate = lambda p: vec3(O[1].x + S.x * (O[0].dot(X) - p.dot(X)),
                                O[1].y + S.y * (O[0].dot(Y) - p.dot(Y)),
                                O[1].z)
        uv = list(locate(self.vertices[v]) for v in self.faces[f])
        uvs[f] = uv
        for k, (i, j) in enumerate(slide(self.faces[f], 2)):
            if not (i, j) in seams:
                adj = self.e2f.get((j, i))
                if adj is not None and adj not in uvs:
                    q = quat.uu(N, face_normals[adj])
                    O = (self.vertices[i], uv[k])
                    self.unwrap_uvs(adj, O, X.cp().rot(q), Y.cp().rot(q), S, seams, uvs)
        return uvs

    def project_uvs_xy(self, s=1.0):
        """Get trivial vertex UVs using xy projection of vertices"""
        S = vec3.U(s)
        uvs = {}
        for f, face in self:
            uvs[f] = [self.vertices[v].xy() * S for v in face]
        return uvs

    def angle_seams(self, alpha=(np.pi / 2)):
        """Find the set of edges joining faces with normals
        differing in angle by alpha or more radians"""
        seams = set()
        fN = self.face_normals()
        for f, face in self:
            for i, j in slide(face, 2):
                #left  = self.e2f.get((i, j))
                #assert left == f # occurs when topology is imperfect
                right = self.e2f.get((j, i))
                if right is None:
                    seams.add((i, j))
                else:
                    if near(fN[f].ang(fN[right]), alpha) >= alpha:
                        seams.add((i, j))
        return seams

    def offset(self, r=1):
        """Contract the surfaces of self along their normals"""
        face_normals = self.face_normals()
        delta = defaultdict(list)
        for v, p in enumerate(self.vertices):
            if p is not None:
                for o in self.v2f[v]:
                    N = face_normals[o]
                    for dN in delta[v]:
                        if dN.isnear(N):
                            break
                    else:
                        delta[v].append((N * -r))
                delta[v] = vec3.sum(delta[v])
                #delta[v] = vec3.sum([face_normals[o] for o in self.v2f[v]]).nrm() * -r
        for v, dp in delta.items():
            self.vertices[v].trn(dp)
        return self

    def dissolve(self):
        """Merge duplicate vertices"""
        groups = defaultdict(list)
        for i, v in enumerate(self.vertices):
            for j in groups:
                if self.vertices[j].isnear(v):
                    groups[i].append(j)
                    print('dupe!')
                    break
            else:
                groups[i].append(i)
        raise NotImplementedError(f'{groups}')

    def _fp(self, p, e=0.00001, min_index=0):
        """Find vertex in neighborhood of p if one exists"""
        for i, o in enumerate(self.vertices):
            if i < min_index:
                continue
            if p.isnear(o, e):
                return i

    def nearest(self, p):
        """Find vertex which is nearest to p"""
        nearest = None
        for i, o in enumerate(self.vertices):
            op = o.d(p)
            if nearest is None or op < nearest[1]:
                nearest = i, op
        if nearest is not None:
            return nearest[0]

    def av(self, p, e=0.00001, min_index=0):
        """Add vertex to the mesh without connectivity"""

        #e = 0.00001 ###

        v = None if e is None else self._fp(p, e, min_index)
        if v is None:
            v = len(self.vertices)
            self.vertices.append(p)
        return v

    def findfaces(self, meta):
        """Find all faces with matching meta data"""
        return [f for f, face in self if (self.meta[f] == meta)]

    def af(self, loop, meta=None, e=0.00001, raise_on_duplicate=True):
        """Add a face to the mesh via a loop of vertices"""
        if not isinstance(loop[0], int):
            loop = [self.av(p, e=e) for p in loop]
        f = len(self.faces)
        for i, j in slide(loop, 2):
            if (i, j) in self.e2f:
                if raise_on_duplicate:
                    raise ValueError
            else:
                self.v2f[i].add(f)
                self.e2f[(i, j)] = f
                self.meta[f] = meta
        self.faces.append(loop)
        self.nface += 1
        return f

    def rf(self, f):
        """Remove face f including its connectivity"""
        for i, j in slide(self.faces[f], 2):
            del self.e2f[(i, j)]
            self.v2f[i].remove(f)
        self.faces[f] = None
        del self.meta[f]
        self.nface -= 1

    def apy(self, py, e=0.00001, h=None, r=10000, **kws):
        """Add triangles covering a polygon, possibly with holes"""
        eloop, iloops = py
        n = loop_normal(eloop)
        q = quat.toxy(n)
        q.rot(eloop)
        for iloop in iloops:
            q.rot(iloop)
        t = triangulation(py, e, h, r)
        q = q.fp()
        q.rot(t.points)
        return [self.af(tri, e=e, **kws) for tri in t.simplices()]

    def fan(self, center, rim, e=0.00001, **kws):
        """Add triangles on surface of cone"""
        new_faces = []
        for u, v in slide(rim, 2, 0):
            new_faces.append(self.af([center, u, v], e=e, **kws))
        return new_faces

    def bridge(self, left, right, m=0, e=0.00001, **kws):
        """Add quadrilaterals between a pair of loops"""
        assert len(left) == len(right)
        if not isinstance(left[0], int):
            left = [self.av(p, e=e) for p in left]
            right = [self.av(p, e=e) for p in right]
        new_faces = []
        for (a, b), (c, d) in zip(slide(left, 2, m), slide(right, 2, m)):

            new_faces.append(self.af([a, b, d, c], e=e, **kws))

            #fan = [a, b, d, c]
            #com = self.av(vec3.com([self.vertices[x] for x in fan]), e=e)
            #new_faces.append(self.fan(com, fan, e=e, **kws))

        return new_faces

    def bridges(self, loops, m=0, n=0, e=0.00001, **kws):
        """Bridge a sequence of loops"""
        if not isinstance(loops[0][0], int):
            loops = [[self.av(p, e=e) for p in loop] for loop in loops]
        new_patches = []
        for u, v in slide(loops, 2, m):
            new_patches.append(self.bridge(u, v, m=n, e=e, **kws))
        return new_patches

    def grid(self, a, b, c, d, n=3, m=3, e=0.00001, **kws):
        """Create grid patch over quadrilateral region"""
        rails = [([d] + d.line(c, n - 1) + [c]), ([a] + a.line(b, n - 1) + [b])]
        lines = [([x] + x.line(y, m - 1) + [y]) for x, y in zip(*rails)]
        return self.bridges(lines, m=1, n=1, e=e, **kws)

    def extrude(self, faces, offset, meta=None):
        """Extrude a set of faces along a vector"""

        # the faces induce a set of seams
        # faces induce cover on 1-1 mappings of modified seams...
        # extrude each edge of each seam using some local/global strategy

        bottom = self.perimeter(faces)
        assert len(bottom) == 1, 'cannot extract surfaces with genus > 0'
        surface = [[self.vertices[v].cp() for v in self.faces[f]] for f in faces]
        for f in faces:
            if meta is None:
                meta = self.meta[f]
            self.rf(f)
        patches = [[self.af(offset.trnps(new), meta=meta) for new in surface]]
        top = self.perimeter(patches[0])
        patches.append(self.bridge(
            [i for i, j in bottom[0]],
            [i for i, j in top[0]], meta=meta))
        return patches

    def revolve(self, curve, axis, n=4):
        patches = []
        back = [p.cp() for p in curve]
        q = quat.av(np.pi * (2 / n), axis)
        for i in range(n):
            front = [p.cp() for p in back]
            q.rot(front)
            patches.append(self.bridge(front, back, m=1))
            back = front
        return patches

    def opposite(self, i, j):
        """Select the opposing edge"""
        f = self.e2f.get((i, j))
        if f is None:
            return
        else:
            face = self.faces[f]
            # is this an implicit quad-only assumption?
            k = face.index(i) - 2
            return face[k], face[k + 1]

    def vertexloop(self, i, j, loop=None):
        """Select a sequence of consecutive edges without turning"""
        raise

    def edgeloop(self, i, j, loop=None):
        """Select a sequence of face-opposing edges"""
        if loop and (i, j) == loop[0]:
            return loop
        else:
            if loop is None:
                loop = []
            loop.append((i, j))
            far = self.opposite(j, i)
            return loop if far is None else self.edgeloop(*far, loop=loop)

    def edgesplit(self, i, j, a=0.5, e=0.00001):
        """Split the face incident to edge ij"""
        new_faces = []
        f = self.e2f.get((i, j))
        if f is not None:
            u = self.av(self.vertices[i].lerp(self.vertices[j], a), e=e)
            face = self.faces[f]
            meta = self.meta[f]
            self.rf(f)
            v = face.index(i) - 2
            if face[v] == j:
                # is this path for tri-only?
                new_faces.append(self.af([i, u, face[v + 1]], meta=meta))
                new_faces.append(self.af([u, j, face[v + 1]], meta=meta))
            else:
                # is this path for quad-only?
                w = face[v + 1]
                v = face[v]
                x = self.av(self.vertices[w].lerp(self.vertices[v], a), e=e)
                new_faces.append(self.af([i, u, x, w], meta=meta))
                new_faces.append(self.af([u, j, v, x], meta=meta))
        return new_faces

    def edgeloopsplit(self, i, j, a=0.5, loop=None):
        """Split the edges incident to each of an edgeloop ij"""
        split = []
        for u, v in self.edgeloop(i, j, loop=loop):
            split.extend(self.edgesplit(u, v, a))
        return split

    def facesplit(self, f, u, v, e=0.00001):
        """Break face in two based on intersection with uv"""
        raise

    def perimeter(self, faces):
        """Find the set of sequences of edges bounding a set of faces"""
        edges = []
        for f in faces:
            for i, j in slide(self.faces[f], 2):
                if (j, i) in edges:
                    edges.remove((j, i))
                else:
                    edges.append((i, j))
        seams = []
        if edges:
            seams = [[edges.pop(0)]]
            while edges:
                for k, (i, j) in enumerate(edges):
                    if i == seams[-1][-1][1]:
                        seams[-1].append((i, j))
                        edges.pop(k)
                        break
                else:
                    edges.pop(k)
                    seams.append([(i, j)])
        return seams
