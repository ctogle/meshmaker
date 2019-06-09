from .vec3 import vec3
from .geometry import sintsxyp
from collections import defaultdict

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


class planargraph:

    def __init__(self, segs=None, epsilon=0.00001):
        self.vertices = []
        self.vertex_count = 0
        self.rings = defaultdict(set)
        self.edges = []
        self.edge_count = 0
        self.edge_lookup = {}
        if segs:
            for u, v in segs:
                self.ae(u, v, epsilon)


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

    def fv(self, p, e=0.00001):
        """find vertex"""
        for i, v in enumerate(self.vertices):
            #if v.isnear(p, e):
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
                if p.insxy(u, v):
                    self.se(i, j, new_vertex)
                    #break
        return new_vertex


    def ne(self, i, j, **properties):
        """new edge"""
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
        self.edges[self.edge_lookup[(i, j)]] = None
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
            return self.edge_lookup[(i, j)]
        else:
            u = self.vertices[i]
            v = self.vertices[j]
            # if u, v intersects existing vertex - ae both segs
            for k, w in enumerate(self.vertices):
                if (not (w.isnear(u, e) or w.isnear(v, e))) and w.insxy(u, v):
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
                ips = sorted(ips, key=lambda ip: ip.dot(u.tov(v)))
                ips = [self.nv(ip) for ip in ips]
                new_edges = []
                for (n, m), ip in zip(ies, ips):
                    new_edges.append(self.se(n, m, ip))
                ips.insert(0, i)
                ips.append(j)
                for k in range(1, len(ips)):
                    new_edges.append(self.ne(ips[k - 1], ips[k], **properties))
                return new_edges
            else:
                return self.ne(i, j, **properties)


    def destem(self):
        stems = []
        for i in self.rings:
            if len(self.rings[i]) == 1:
                stems.append(i)
        if stems:
            for i in stems:
                self.rv(i)
            return self.destem()

    def minlen(self, epsilon=None):
        """"""
        lengths = []
        for i in self.rings:
            for j in self.rings[i]:
                lengths.append(self.vertices[i].tov(self.vertices[j]).mag())
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


    def follow(self, i, j, z):
        ring = tuple(self.rings[j])
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
        loop.pop()
        lowest = loop.index(min(loop))
        return tuple(loop[lowest:] + loop[:lowest])

    def loops(self):
        loops = set()
        for i, j, properties in filter(lambda e: bool(e), self.edges):
            loops.add(self.loop(i, j, -1))
        return list(loops)


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


