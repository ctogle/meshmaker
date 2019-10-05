#from dilap.geometry import vec3, pointset
#from dilap.geometry.tools import orient2d, incircle, circumscribe_tri
#from dilap.geometry.tools import orient2d, incircle
#from dilap.geometry.polymath import sintsxy, sintsxyp, tinbxy, circumscribe_tri

from .plt import plt, plot, plot_loop, plot_point
from .vec3 import vec3
from .geometry import orient2d
from .geometry import sintsxyp, loop_contains_triangle, circumscribe_triangle

import numpy as np
from tqdm import tqdm


class triangulation:

    def simplices(self):
        """yield each triangle"""
        for j, t in enumerate(self.triangles):
            if t:
                yield tuple(self.points[x].cp() for x in t)

    def onboundary(self, v, e):
        """determine if vertex v is on the boundary (ghost edges)"""
        p = self.points[v]
        for j, ghost in enumerate(self.ghosts):
            if not ghost is None:
                u, v, _ = ghost
                up, vp = self.points[u], self.points[v]
                if p.onsxy(up, vp, 1):
                    return j
        else:
            return -1

    def adjacent(self, u, v):
        """find vertex w such that uv is a positively oriented edge"""
        if (u, v) in self.trianglelookup:
            triangle = self.trianglelookup[(u, v)]
            if triangle is not None:
                x, y, z = self.triangles[triangle]
                if (u, v) == (y, z) or (u, v) == (z, y):
                    return x
                elif (u, v) == (z, x) or (u, v) == (x, z):
                    return y
                elif (u, v) == (x, y) or (u, v) == (y, x):
                    return z
                else:
                    raise ValueError
        if (u, v) in self.ghostlookup:
            ghost = self.ghostlookup[(u, v)]
            return -1 if ghost is None else -2
        return -1

    def locallydelaunay(self, u, v):
        """determine if uv is a delaunay edge"""
        x, y = self.adjacent(u, v), self.adjacent(v, u)
        if not any((x == -1, x == -2, y == -1, y == -2)):
            up, vp = self.points[u], self.points[v]
            xp, yp = self.points[x], self.points[y]
            #if sintsxy(up, vp, xp, yp, col=0):
            ip = sintsxyp(up, vp, xp, yp, colinear=0)
            if ip:
                #if incircle(up, vp, xp, yp) > 0:
                if yp.incircle(up, vp, xp) > 0:
                    return False
                #if incircle(vp, up, yp, xp) > 0:
                if xp.incircle(vp, up, yp) > 0:
                    return False
        return True

    def plot(self, ax=None, hl=[]):
        #ax = ax if ax else plot_axes_xy(100)
        plot_loop(ax, self.polygon[0], lw=2, col='m')
        for hole in self.polygon[1]:
            plot_loop(ax, hole, lw=2, col='y')
        for triangle in self.triangles:
            if triangle:
                #u, v, w = self.points.gps(triangle)
                u, v, w = tuple(self.points[x].cp() for x in triangle)
                #cc, cr = circumscribe_tri(u, v, w)
                #plot_circle_xy(cc, cr, ax, lw=2, col='k')
                plot_loop(ax, [u, v, w], lw=1, col='g')
        for i in hl:
            triangle = self.triangles[i]
            if triangle:
                u, v, w = tuple(self.points[x].cp() for x in triangle)
                plot_loop(ax, [u, v, w], lw=2, col='y')
        #for ghost in self.ghosts:
        #    if ghost:
        #        plot_edges_xy(self.points.gps((ghost[0], ghost[1])), ax, lw=3, col='b')
        for i, p in enumerate(self.points):
            if i > 3:
                plot_point(ax, p, col='r')
        return ax

    def __init__(self, polygon, epsilon, hmin=2, crad=10000):
        #self.points = pointset()
        self.points = []
        self.ntriangles, self.nghosts = 0, 0
        self.triangles, self.ghosts = [], []
        self.trianglelookup, self.ghostlookup = {}, {}

        self.epsilon = epsilon
        self.polygon = polygon
        splitpolygon = triangulation.delaunayedges(self.polygon, self.epsilon)
        if hmin:
            hmin, splitpolygon = triangulation.chewfirstedges(
                             splitpolygon, self.epsilon, hmin)
        self.hmin = hmin
        self.splitpolygon = splitpolygon
        self.coverradius = crad

        self.cover = self.initialize(self.splitpolygon[0], self.epsilon)
        self.edges = []
        self.edges.extend(self.insertloop(self.splitpolygon[0], self.epsilon))
        for hole in self.splitpolygon[1]:
            self.edges.extend(self.insertloop(hole, self.epsilon))
        self.prune(self.splitpolygon[0], self.splitpolygon[1], self.edges, self.epsilon)
        self.delaunay()
        if self.hmin:
            print('refine with hmin: %0.3f' % self.hmin)
            self.refine(self.hmin, self.epsilon)
        # smooth laplacian

    @staticmethod
    def delaunayedges(py, e):
        """return copy of a polygon with potentially additional
        points such that each edge is locally delaunay"""
        points = py[0][:] + [x for y in py[1] for x in y]

        def splitloop(loop):
            for j in range(len(loop)):
                u, v = loop[j - 1], loop[j]
                w = u.lerp(v, 0.5)
                r = u.d(w)
                for p in points:
                    #print(u.d(p), v.d(p))
                    if not (u.d(p) < e or v.d(p) < e):
                        if w.d(p) < r:
                            loop.insert(j, w)
                            points.append(w)

                            #print(loop)
                            #print('-' * 50)
                            #f, ax = plot()
                            #plot_loop(ax, loop, lw=1, mk='o')
                            #plot_loop(ax, w.ring(r, 16))
                            #plt.show()

                            return True

        def handleloop(loop):
            done = False
            while not done:
                done = not splitloop(loop)
            return loop

        return handleloop(py[0][:]), [handleloop(h[:]) for h in py[1]]

    @staticmethod
    def chewfirstedges(py, e, hmax):
        """split edges to prepare for chews first refinement algorithm"""
        eb, ibs = py
        els = [eb[x-1].d(eb[x]) for x in range(len(eb))]
        for ib in ibs:
            els.extend([ib[x-1].d(ib[x]) for x in range(len(ib))])

        hmin = min(els) * np.sqrt(3)
        hmin = min(hmax, hmin)

        def splitloop(loop):
            o = []
            l = loop[-1]
            for j in range(len(loop)):
                n = loop[j]
                el = n.d(l)
                m = 1
                while el / m > hmin:
                    m += 1
                ps = l.line(n, m - 1)
                if ps:
                    for dp in ps:
                        o.append(dp)
                o.append(n)
                l = n
            return o

        splitpy = (splitloop(eb), [splitloop(h) for h in ibs])
        return hmin, splitpy

    def _fp(self, p, e):
        """search for known points"""
        for i, q in enumerate(self.points):
            if q.isnear(p, e=e):
                return i
        else:
            return -1

    def _ap(self, p):
        """add a new point"""
        i = len(self.points)
        self.points.append(p)
        return i

    def _nps(self, *ps, e=0.00001):
        """add potentially new points to the point set"""
        new = []
        for np in ps:
            i = self._fp(np, e=e)
            i = self._ap(np) if i == -1 else i
            new.append(i)
        return new

    def initialize(self, boundary, e):
        """provide an initial triangle cover for the boundary"""
        u, v, w, x = vec3.com(boundary).ring(self.coverradius, 4)
        u, v, w, x = self._nps(u, v, w, x, e=e)
        self.addtriangle(u, v, w)
        self.addtriangle(u, w, x)
        self.addghost(v, u)
        self.addghost(w, v)
        #self.addghost(u, w)
        self.addghost(x, w)
        self.addghost(u, x)
        return (u, v, w, x)

    def insertloop(self, loop, e):
        """perform point location for each point in a loop of points"""
        npoints = len(loop)
        edges = [(loop[j - 1], loop[j]) for j in range(npoints)]
        #loop = sorted(loop)
        for j in tqdm(range(npoints), total=npoints):
            self.insertpoint(loop[j], e)
        return edges

    def insertpoint(self, p, e):
        """point location modifies triangulation to include new point p"""
        oc = self.ntriangles
        nv = self._fp(p, e=e)
        if nv == -1:
            nv = self._ap(p)
            g = self.onboundary(nv, e)
            if g == -1:
                for j, t in enumerate(self.triangles):
                    if t:
                        u, v, w = tuple(self.points[x] for x in t)
                        #a, b, c = orient2d(p, u, v), orient2d(p, v, w), orient2d(p, w, u)
                        #if triangulation.pintxy(p, u, v, w, e):
                        #if not any((a < 0, b < 0, c < 0)):
                        #if p.intrixy(u, v, w, e):
                        if p.inbxy((u, v, w)):
                            self.insertvertex(nv, *t)
                            return range(oc, self.ntriangles)
                else:
                    for j, g in enumerate(self.ghosts):
                        if not g is None:
                            u, v, w = g
                            #up, vp = self.points.gps((u, v))
                            up, vp = self.points[u], self.points[v]
                            if orient2d(up, vp, p) >= 0:
                                self.insertghostvertex(nv, u, v, w, e)
                                #raise ValueError
                                return range(oc, self.ntriangles)
            else:
                u, v, _ = self.ghosts[g]
                self.deleteghost(u, v)
                self.addghost(u, nv)
                self.addghost(nv, v)
                w = self.adjacent(v, u)
                self.deletetriangle(v, u, w)
                self.addtriangle(w, nv, u)
                self.addtriangle(w, v, nv)
                raise ValueError
                return range(oc, self.ntriangles)
            raise ValueError('failed to locate point!')
        else:
            print(self.points[nv], p, self.points[nv].d(p))
            raise ValueError('point already located!')

    def prune(self, boundary, holes, edges, e):
        '''
        remove triangles either outside of the boundary or inside one of the holes
        replace ghost triangles such that the targeted edges are included
        '''
        extras = []
        for triangle in self.triangles:
            if triangle:
                #up, vp, wp = self.points.gps(triangle)
                up, vp, wp = tuple(self.points[x] for x in triangle)
                #if tinbxy(up, vp, wp, boundary):
                if loop_contains_triangle(boundary, up, vp, wp):
                    for hole in holes:
                        #if tinbxy(up, vp, wp, hole):
                        if loop_contains_triangle(hole, up, vp, wp):
                            extras.append(triangle)
                            break
                else:
                    extras.append(triangle)
        for triangle in extras:
            self.deletetriangle(*triangle)
        for j, g in enumerate(self.ghosts):
            if g:
                u, v, _ = g
                self.deleteghost(u, v)
        for u, v in edges:
            #u, v = self.points.fps(edge, e)
            u, v = self._fp(u, e=e), self._fp(v, e=e)
            if self.adjacent(u, v) == -1:
                self.addghost(u, v)
            if self.adjacent(v, u) == -1:
                self.addghost(v, u)

    def insertvertex(self, u, v, w, x):
        '''
        insert new vertex u which lies inside of triangle vwx
        '''
        #up, vp, wp, xp = self.points.gps((u, v, w, x))
        #up, vp = self.points[u], self.points[v]
        #wp, xp = self.points[w], self.points[x]
        self.deletetriangle(v, w, x)
        self.digcavity(u, v, w)
        self.digcavity(u, w, x)
        self.digcavity(u, x, v)

    def insertghostvertex(self, u, v, w, x, e):
        '''
        insert new ghost vertex u with bordering edge vw
        '''
        assert(x == -2)
        self.deleteghost(v, w)
        self.addghost(v, u)
        self.addghost(u, w)
        if not self.onboundary(u, e) == -1:
            self.addtriangle(u, v, w)

    def addtriangle(self, u, v, w):
        '''
        insert new triangle uvw
        '''
        self.triangles.append((u, v, w))
        self.trianglelookup[(u, v)] = self.ntriangles
        self.trianglelookup[(v, w)] = self.ntriangles
        self.trianglelookup[(w, u)] = self.ntriangles
        self.ntriangles += 1

    def deletetriangle(self, u, v, w):
        '''
        remove triangle uvw
        '''
        which = self.trianglelookup[(u, v)]
        if which is not None:
            self.triangles[which] = None
        self.trianglelookup[(u, v)] = None
        self.trianglelookup[(v, w)] = None
        self.trianglelookup[(w, u)] = None

    def addghost(self, u, v):
        '''
        insert new ghost triangle adjacent to uv
        '''
        self.ghosts.append((u, v, -2))
        self.ghostlookup[(u, v)] = self.nghosts
        self.nghosts += 1

    def deleteghost(self, u, v):
        '''
        remove ghost triangle adjacent to uv
        '''
        which = self.ghostlookup[(u, v)]
        if which is not None:
            self.ghosts[which] = None
        self.ghostlookup[(u, v)] = None

    def digcavity(self, u, v, w):
        '''
        determine if potential new triangle uvw is delaunay
        u is a new vertex; vw is an existing edge
        '''
        up, vp, wp = self.points[u], self.points[v], self.points[w]
        x = self.adjacent(w, v)
        if x == -1:
            return
        elif x == -2:
            # adjacent triangle is a ghost
            self.addtriangle(u, v, w)
        else:
            xp = self.points[x]

            '''
            ax = plot_axes_xy(125, (0, 0))
            ax = self.plot(ax)
            for j, p in enumerate(self.points.ps):
                plot_point_xy(p, plot_point_xy_annotate(p, ax, str(j)), col='g', mk='*')
            plot_polygon_xy((up, vp, wp), ax, col='r', lw=2)
            plot_point_xy(up, plot_point_xy_annotate(up, ax, 'up    '), col='y', mk='s')
            plot_point_xy(vp, plot_point_xy_annotate(vp, ax, 'vp    '), col='m', mk='*')
            plot_point_xy(wp, plot_point_xy_annotate(wp, ax, 'wp    '), col='m', mk='*')
            plot_point_xy(xp, plot_point_xy_annotate(xp, ax, 'xp    '), col='m', mk='*')
            cc, cr = circumscribe_tri(up, vp, wp)
            plot_circle_xy(cc, cr, ax, lw=4, col='r')

            plt.show()
            print('dug cavity', 'incircle', incircle(up, vp, wp, xp))
            '''

            #if incircle(up, vp, wp, xp) > 0:
            if xp.incircle(up, vp, wp) > 0:
                # u, v, w is not delaunay; dig adjacent triangles
                self.deletetriangle(w, v, x)
                self.digcavity(u, v, x)
                self.digcavity(u, x, w)
            else:
                # u, v, w is delaunay
                self.addtriangle(u, v, w)

    def flipedge(self, u, v):
        '''
        flip the edge uv (replace 2 triangles with two different triangles)
        uv must bound two non-ghost triangles
        '''
        x, y = self.adjacent(u, v), self.adjacent(v, u)
        self.deletetriangle(u, v, x)
        self.deletetriangle(v, u, y)
        self.addtriangle(x, y, v)
        self.addtriangle(y, x, u)
        return (y, v), (v, x), (x, u), (u, y)

    def delaunay(self):
        '''
        ensure triangulation is delaunay via edge flips
        '''
        unfinished = list(self.trianglelookup.keys())
        while unfinished:
            u, v = unfinished.pop(0)
            if not self.locallydelaunay(u, v):
                newedges = self.flipedge(u, v)
                unfinished.extend(newedges)

    def refine(self, h, e):
        '''
        perform chew's first refinement algorithm
        '''


        '''
        print('prerefinement!')
        ax = plot_axes_xy(10, (-50, 50))
        self.plot(ax)
        plt.show()
        '''
        cover = tuple(self.points[x] for x in self.cover)
        unfinished = [t for t in self.triangles if t is not None]
        while unfinished:
            t = unfinished.pop(0)
            if not t in self.triangles:
                continue
            #u, v, w = self.points.gps(t)
            u, v, w = tuple(self.points[x] for x in t)
            #tcp, tcr = circumscribe_triangle(u, v, w, e)
            tcp, tcr = circumscribe_triangle(u, v, w)
            if tcr > h:

                #oob = not tcp.inbxy(self.points.gps(self.cover))
                oob = not tcp.inbxy(cover)
                oob = oob or not tcp.inbxy(self.polygon[0])
                oob = oob or any((tcp.inbxy(h) for h in self.polygon[1]))
                if oob:
                    print('OOOB', oob)
                    '''
                    print('refinement insertion!', tcr, tcr / h)
                    ax = plot_axes_xy(125, (0, 0))
                    self.plot(ax)
                    #for j, t in enumerate(self.triangles):
                    #    if t:
                    #        t = self.points.gps(t)
                    #        plot_polygon_xy(t, ax, col='g', lw=4)
                    plot_polygon_xy((u, v, w), ax, col='m', lw=6)
                    plot_circle_xy(tcp, tcr, ax, center=True, col='k', lw=4)
                    plt.show()

                    break
                    '''

                newtriangles = self.insertpoint(tcp, e)
                unfinished.extend([self.triangles[j] for j in newtriangles])

                #print('refinement insertion!', tcr, tcr / h)
                #ax = plot_axes_xy(20, (-48, 55))
                #self.plot(ax)
                ##for j, t in enumerate(self.triangles):
                ##    if t:
                ##        t = self.points.gps(t)
                ##        plot_polygon_xy(t, ax, col='g', lw=4)
                #plot_polygon_xy((u, v, w), ax, col='m', lw=6)
                #plot_circle_xy(tcp, tcr, ax, center=True, col='k', lw=4)
                #plt.show()

