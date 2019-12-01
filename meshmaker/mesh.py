import numpy as np
from collections import defaultdict
from .vec3 import vec3
from .quat import quat
from .delaunay import triangulation
from .geometry import sintsxyp, loop_area, loop_normal, slide, isnear

from .plt import plot, plot_pg, plot_point, plot_edge


class trimesh:

    def __init__(self):
        self.vertices = []
        self.faces = []

    def __iter__(self):
        for i, j, k in (f for f in self.faces if f):
            u = self.vertices[i][0]
            v = self.vertices[j][0]
            w = self.vertices[k][0]
            yield u, v, w

    def av(self, v):
        self.vertices.append(v)
        return len(self.vertices) - 1

    def af(self, p1, p2, p3, p4=None):
        nm = (p1.tov(p2)).crs(p1.tov(p3)).nrm()
        u1 = vec3(0, 0, 0)
        u2 = vec3(1, 0, 0)
        u3 = vec3(1, 1, 0)
        v1 = self.av((p1, nm.cp(), u1))
        v2 = self.av((p2, nm.cp(), u2))
        v3 = self.av((p3, nm.cp(), u3))
        if p4:
            u4 = vec3(0, 1, 0)
            v4 = self.av((p4, nm.cp(), u4))
            self.faces.append((v1, v2, v3))
            self.faces.append((v1, v3, v4))
        else:
            self.faces.append((v1, v2, v3))

    def apy(self, py, e=0.00001, h=None, r=10000):
        eloop, iloops = py
        n = loop_normal(eloop)
        q = quat.toxy(n)
        q.rot(eloop)
        for iloop in iloops:
            q.rot(iloop)
        t = triangulation(py, e, h, r)
        q = q.fp()
        q.rot(t.points)
        for x, y, z in t.simplices():
            self.af(x, y, z)

class planargraph:


    def merge(self, o):
        vmap = {}
        for i, v in enumerate(o.vertices):
            if v is None:
                vmap[i] = len(self.vertices)
                self.vertices.append(None)
            else:
                vmap[i] = self.nv(v.cp(), **v.properties)
        emap = {}
        for i, e in enumerate(o.edges):
            if e is None:
                emap[i] = len(self.edges)
                self.edges.append(None)
            else:
                u, v, properties = e
                emap[i] = self.ne(vmap[u], vmap[v], **properties)
        return vmap, emap


    def cp(self):
        o = self.__class__()
        for v in self.vertices:
            if v is None:
                o.vertices.append(None)
            else:
                o.nv(v.cp(), **v.properties)
        for e in self.edges:
            if e is None:
                o.edges.append(None)
            else:
                i, j, e = e
                o.ne(i, j, **e)
        return o


    def __init__(self, segs=None, epsilon=0.00001):
        self.vertices = []
        self.vertex_count = 0
        self.rings = defaultdict(set)
        self.edges = []
        self.edge_count = 0
        self.edge_lookup = {}
        if segs:
            for u, v in segs:
            #    i, j = self.av(u), self.av(v)
            #    self.ne(i, j)
            #self.dedupe()
                try:
                    self.ae(u, v, epsilon)
                except AssertionError:
                    #print('segment too short')
                    pass


    def nv(self, p, **properties):
        """new vertex"""
        p.properties = properties
        new_vertex = len(self.vertices)
        self.vertices.append(p)
        self.vertex_count += 1
        return new_vertex

    def rv(self, i):
        """remove vertex"""
        self.vertices[i] = None
        self.vertex_count -= 1
        for j in self.rings[i].copy():
            self.re(i, j)
        del self.rings[i]

    def sv(self, i, alpha):
        """split vertex - insert new vertex at location of vertex i,
        attach edge to existing vertex i and reconnect all edges to
        the left of the plane defined by x-axis signed angle alpha"""
        plane = vec3(np.cos(alpha), np.sin(alpha), 0)
        j = self.nv(self.vertices[i].cp())
        for k in tuple(self.rings[i]):
            uv = (self.vertices[k] - self.vertices[i]).nrm()
            if plane.saxy(uv) < np.pi:
                self.re(i, k)
                self.ne(k, j)
        self.ne(i, j)
        return j

    def fv(self, p, e=0.00001):
        """find vertex"""
        for i, v in enumerate(self.vertices):
            #if v.isnear(p, e):
            if v is not None:
                if v.dxy(p) < e:
                    return i

    def av(self, p, e=0.00001, **properties):
        """add vertex

        if p is near an existing vertex, return that vertex
        if p is near an existing edge, bisect that edge and return new vertex
        """
        new_vertex = self.fv(p, e) if e > 0 else None
        if new_vertex is None:
            new_vertex = self.nv(p, **properties)
            for i, j in tuple(self.edge_lookup.keys()):
                u, v = self.vertices[i], self.vertices[j]
                if p.insxy(u, v, e=e):
                    self.se(i, j, new_vertex)
                    #break
        return new_vertex


    def ne(self, i, j, allow_loop=False, **properties):
        """new edge"""
        if i == j and not allow_loop:
            return
        i, j = (i, j) if (i < j) else (j, i)
        self.rings[i].add(j)
        self.rings[j].add(i)
        new_edge = len(self.edges)
        self.edges.append((i, j, properties))
        self.edge_lookup[(i, j)] = new_edge
        self.edge_count += 1
        return new_edge

    def re(self, i, j):
        """remove edge"""
        i, j = (i, j) if (i < j) else (j, i)
        k = self.edge_lookup.get((i, j))
        if k is not None:
            self.edges[k] = None
            del self.edge_lookup[(i, j)]
            self.rings[i].remove(j)
            self.rings[j].remove(i)
            self.edge_count -= 1

    def se(self, i, j, k):
        """split edge"""
        i, j = (i, j) if (i < j) else (j, i)
        _, _, properties = self.edges[self.edge_lookup[(i, j)]]
        self.re(i, j)
        u = self.ne(i, k, **properties)
        v = self.ne(k, j, **properties)
        return u, v

    def ae(self, i, j, e=0.00001, **properties):
        """add edge

        if (i, j) is an existing edge, return that edge
        if (i, j) intersects an existing vertex, add split edges
        if (i, j) intersects another edge skew and without any endpoints
            split intersecting edges and 4 new edges

        """
        assert i != j
        i = i if isinstance(i, int) else self.av(i, e)
        j = j if isinstance(j, int) else self.av(j, e)
        i, j = (i, j) if (i < j) else (j, i)
        assert i != j
        if (i, j) in self.edge_lookup:
            # edge already exists; return its index
            return self.edge_lookup[(i, j)]
        else:
            # find intersections with existing geometry
            u = self.vertices[i]
            v = self.vertices[j]

            # if u, v intersects existing vertex - ae both segs
            for k, w in enumerate(self.vertices):
                if w is not None:
                    if (not (w.isnear(u, e) or w.isnear(v, e))) and w.insxy(u, v, e):
                        return self.ae(i, k, e), self.ae(k, j, e)
            # 
            ips, ies = [], []
            for n, m in self.edge_lookup:
                p, q = self.vertices[n], self.vertices[m]
                ip = sintsxyp(u, v, q, p, False, False, False, True)
                if ip:
                    ips.append(ip)
                    ies.append((n, m))

            if ips:
                # add new edges forming path between ips; return indices
                zipped = zip(ips, ies)
                zipped = sorted(zipped, key=lambda ip: ip[0].dot(u.tov(v)))
                new_edges = []
                for k, (ip, (n, m)) in enumerate(zipped):
                    ip = self.nv(ip)
                    new_edges.append(self.se(n, m, ip))
                    ips[k] = ip
                ips.insert(0, i)
                ips.append(j)
                for k in range(1, len(ips)):
                    new_edges.append(self.ne(ips[k - 1], ips[k], **properties))
                return new_edges
            else:
                # add new edge; return its index
                return self.ne(i, j, **properties)

    def dedupe(self):
        """deduplicate vertices and edges (adds intersections)"""
        unique = set()
        for i, p in enumerate(self.vertices):
            if p is not None:
                for j in unique:
                    q = self.vertices[j]
                    if p.isnear(q):
                        for k in self.rings[i]:
                            if not j == k:
                                self.ne(j, k)
                        self.rv(i)
                        break
                else:
                    unique.add(i)
        for i, u in enumerate(self.edges):
            if u is not None:
                k, l, _ = u
                p, q = self.vertices[k], self.vertices[l]
                for j, v in enumerate(self.edges):
                    if j >= i:
                        break
                    if v is not None:
                        n, m, _ = v
                        r, s = self.vertices[n], self.vertices[m]
                        ips = sintsxyp(p, q, r, s,
                                       endpoint=False, endpoints=False,
                                       skew=True, colinear=True)
                        if isinstance(ips, (tuple, list)):
                            self.re(k, l)
                            self.ae(k, l)
                            self.re(n, m)
                            self.ae(n, m)
                        elif isinstance(ips, vec3):
                            self.av(ips)

    def destem(self):
        stems = []
        for i in self.rings:
            if len(self.rings[i]) == 1:
                stems.append(i)
        if stems:
            for i in stems:
                self.rv(i)
            return self.destem()

    def dissolve(self):
        for i, v in enumerate(self.vertices):
            if v is not None:
                if len(self.rings[i]) == 2:
                    j, k = tuple(self.rings[i])
                    x, y, z = self.vertices[j], self.vertices[i], self.vertices[k]
                    ccw = (y - x).crs(z - y).z
                    if isnear(ccw, 0):
                        self.rv(i)
                        self.ae(j, k)

    def minlen(self, epsilon=None):
        """"""
        lengths = []
        for i in self.rings:
            for j in self.rings[i]:
                l = (self.vertices[j] - self.vertices[i]).mag()
                lengths.append(l)
        minlen = min(lengths)
        if epsilon and minlen < epsilon:
            raise ValueError(f'bad minlen: {minlen} (epsilon: {epsilon})')
        return minlen

    def forked(self, i):
        """"""
        return not len(self.rings[i]) == 2

    def forkorder(self, i):
        """compute the topological distance to the nearest vertex of valence > 2"""
        valence = len(self.rings[i])
        if valence == 0:
            min_order = 1
        elif valence == 1:
            a, = tuple(self.rings[i])
            min_order = len(self.loop_until(i, a, 1, self.forked)) - 1
        elif valence == 2:
            a, b = tuple(self.rings[i])
            n_a, n_b = len(self.rings[a]), len(self.rings[b])
            if (n_a == 1 and n_b == 1) or (n_a > 2 or n_b > 2):
                min_order = 1
            elif n_a == 1:
                min_order = len(self.loop_until(i, b, 1, self.forked)) - 1
            elif n_b == 1:
                min_order = len(self.loop_until(i, a, 1, self.forked)) - 1
            else:
                a = self.loop_until(i, a, 1, self.forked)
                b = self.loop_until(i, b, 1, self.forked)
                a = len(a) - 1 if len(self.rings[a[-1]]) > 1 else 1e10
                b = len(b) - 1 if len(self.rings[b[-1]]) > 1 else 1e10
                min_order = min(a, b)
        else:
            min_order = 0
        return min_order

    def tangent(self, i, j):
        """Return contiguous, colinear edges to edge (i, j)"""
        i, j = ((j, i) if (j < i) else (i, j))
        wall = [self.edge_lookup[(i, j)]]
        a, b = i, j
        newtip = self.follow_tangential(a, b)
        while newtip is not None:
            key = ((b, newtip) if (b < newtip) else (newtip, b))
            wall.append(self.edge_lookup[key])
            a = b
            b = newtip
            newtip = self.follow_tangential(a, b)
        # TODO: deduplicate redundant code below
        a, b = j, i
        newtip = self.follow_tangential(a, b)
        while newtip is not None:
            key = ((b, newtip) if (b < newtip) else (newtip, b))
            #wall.append(self.edge_lookup[key])
            wall.insert(0, self.edge_lookup[key])
            a = b
            b = newtip
            newtip = self.follow_tangential(a, b)
        return wall

    def follow_tangential(self, i, j):
        """Like `follow` but only moves in a straight line"""
        u = self.vertices[i]
        for k in self.rings[j]:
            if k == i:
                continue
            v = self.vertices[j]
            w = self.vertices[k]
            a = (v - u).saxy((w - v))
            if isnear(a, 0.0, 0.001) or isnear(a, 2 * np.pi, 0.001):
                return k

    def follow(self, i, j, z):
        ring = tuple(self.rings[j])
        assert i in ring, 'must follow an exist edge'
        n_ring = len(ring)
        if n_ring == 1:
            return i
        elif n_ring == 2:
            u, v = ring
            return u if v == i else v
        else:
            alphas = []
            u = self.vertices[i]
            v = self.vertices[j]
            for k in ring:
                if not k == i:
                    w = self.vertices[k]
                    vu, vw = v.tov(u), v.tov(w)
                    alphas.append((vu.saxy(vw), k))
            alphas = sorted(alphas, key=(lambda a: a[0]))
            return alphas[0 if z > 0 else -1][1]

    def loop_until(self, i_0, j_0, z, f):
        loop = [i_0, j_0]
        i, j = i_0, j_0
        while True:
            k = self.follow(i, j, z)
            if f(k):
                loop.append(k)
                break
            elif k == j_0 and j == i_0:
                break
            else:
                loop.append(k)
                i = j
                j = k
        return tuple(loop)

    def loop(self, i_0, j_0, z):
        # TODO: can this use self.loop_until?
        loop = [i_0, j_0]
        i, j = i_0, j_0
        while True:
            k = self.follow(i, j, z)
            if k == j_0 and j == i_0:
                break
            else:
                loop.append(k)
                i = j
                j = k
            if len(loop) > 10000:
                raise ValueError('bad geometry; failed to loop')
        loop.pop()
        lowest = loop.index(min(loop))
        return tuple(loop[lowest:] + loop[:lowest])

    def loops(self):
        loops = set()
        for i, j, properties in filter(lambda e: bool(e), self.edges):
            loops.add(self.loop(i, j, -1))
        return list(loops)

    def polygon(self):
        loops = [[self.vertices[i].cp() for i in l] for l in self.loops()]
        for loop in loops:
            if loop_normal(loop).z < 0:
                loop.reverse()
            assert (not loop_normal(loop).z < 0)
        loops = sorted(loops, key=lambda l: abs(loop_area(l)), reverse=True)
        eloop = loops.pop(0)
        return eloop, loops

    def loop_graph(self):
        """generate a graph where each vertex is associated with a loop,
        each edge with one or more shared loop edges (i.e. adjacent loops),
        and each vertex associates with exactly one loop in loops via an index"""
        pg = planargraph()
        pg.source = self
        loops = self.loops()
        loopps = []
        for i, loop in enumerate(loops):
            ps = [self.vertices[j].cp() for j in loop]
            if loop_normal(ps).z < 0:
                ps.reverse()
            loopps.append((i, loop, vec3.com(ps), loop_area(ps)))
        loopps.sort(key=(lambda v: v[-1]))
        eloop = loopps.pop(-1)
        bound = set()
        for p, q in slide(list(loops[eloop[0]]), 2):
            bound.add((p, q))
            bound.add((q, p))
        loopvs = []
        for i, loop, com, area in loopps:
            nv = pg.nv(com, index=i, loops=[loop], area=area, annotation=' ')
            loopvs.append(nv)
        for i, u in enumerate(loopvs):
            uloop = loops[pg.vertices[u].properties['index']]
            uedges = set()
            for p, q in slide(list(uloop), 2):
                uedges.add((p, q))
                uedges.add((q, p))
            boundary_edges = uedges.intersection(bound)
            pg.vertices[u].properties['boundary'] = boundary_edges
            for j, v in enumerate(loopvs):
                if not i == j:
                    vloop = loops[pg.vertices[v].properties['index']]
                    vedges = set()
                    for p, q in slide(list(vloop), 2):
                        vedges.add((p, q))
                        vedges.add((q, p))
                    shared = uedges.intersection(vedges)
                    if shared:
                        pg.ne(i, j, seam=shared)
        return pg


if __name__ == '__main__':
    pg = planargraph()
    pg.ae(vec3(-1, 0, 0), vec3(1, 0, 0))
    pg.ae(vec3(0, -1, 0), vec3(0, 1, 0))
    print(pg.vertex_count)
    quit()

    pg = planargraph()
    pg.ae(vec3(0, 0, 0), vec3(1, 0, 0))
    pg.ae(vec3(1, 0, 0), vec3(1, 1, 0))
    pg.ae(vec3(1, 1, 0), vec3(0, 1, 0))
    pg.ae(vec3(0, 1, 0), vec3(0, 0, 0))
    pg.ae(vec3(1, 1, 0), vec3(0, 0, 0))
    pg.ae(vec3(1, 1, 0), vec3(2, 2, 0))

    for l in pg.loops():
        print('l', l)

    from .plt import plot, plot_pg
    f, ax = plot()
    plot_pg(ax, pg)
    ax.set_title('%d vertices, %d edges' % (pg.vertex_count, pg.edge_count))
    plt.show()


