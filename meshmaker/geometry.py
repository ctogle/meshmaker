import numpy as np


def near(a, b, epsilon=0.00001):
    d = a - b
    return b if (d * d) < epsilon else a


def isnear(a, b, epsilon=0.00001):
    d = a - b
    return 1 if (d * d) < epsilon else 0


def periodic(a, b, c):
    # TODO: handle c being > 1 period outside of [a, b)
    if c < a:
        return c + (b - a)
    elif c >= b:
        return c - (b - a)
    else:
        return c


def inrng(u, v, w):
    u = near(near(u, v), w)
    return v < u and u < w


def orient2d(a, b, c, epsilon=0.00001):
    """Signed area of the triangle a - c, b - c
        0 if a, b, and c are colinear
    """
    m11 = a.x - c.x
    m12 = a.y - c.y
    m21 = b.x - c.x
    m22 = b.y - c.y
    det = m11 * m22 - m12 * m21
    return near(det, 0, epsilon=epsilon)


def lintp(u, v, p, n):
    """find intersection of line uv and plane p,n"""
    t = (v - u).nrm()
    d = (p - u).dot(n) / t.dot(n)
    return u + t * d


def sintsxyp(u, v, p, q, endpoint=True, endpoints=True, colinear=True, skew=True):
    utov, ptoq, utop = u.tov(v), p.tov(q), u.tov(p)
    uvcrspq, upcrsuv = near(utov.crs(ptoq).z, 0), near(utop.crs(utov).z, 0)
    assert(isnear(utov.z, 0) and isnear(ptoq.z, 0))
    if uvcrspq == 0 and upcrsuv == 0:
        if colinear:
            uv, pq = utov.mag(), ptoq.mag()
            #uvdotpq = uv.dot(pq), up.dot(uv)
            uvdotpq = utov.dot(ptoq)
            updotuv = utop.dot(utov)
            # wtf is with this???
            #assert(not (uvdotpq == 0))
            # wtf is with this???
            t0 = near(near(updotuv / (uv ** 2), 0), 1)
            t1 = near(near(t0 + uvdotpq / (uv ** 2), 0), 1)
            if uvdotpq < 0:
                t0, t1 = t1, t0
            if inrng(t0, 0, 1) or inrng(t1, 0, 1):
                # if include uncontained overlaps
                return sorted((u, v, p, q), key=lambda p: p.dot(utov))[1:-1]
            elif inrng(0, t0, t1) and inrng(1, t0, t1):
                # if include contained overlaps with neither endpoint
                return sorted((u, v, p, q), key=lambda p: p.dot(utov))[1:-1]
            elif (t0 == 0 and t1 == 1) or (t0 == 1 and t1 == 0):
                # if include perfect overlaps
                if endpoints:
                    return u.cp(), v.cp()
            elif ((t0 == 0 or t1 == 1) and t0 < t1) or ((t0 == 1 or t1 == 0) and t0 > t1):
                # if include contained overlaps with one endpoint
                return sorted((u, v, p, q), key=lambda p: p.dot(utov))[1:-1]
                #if endpoint:
                #    return sorted((u, v, p, q), key=lambda p: p.dot(utov))[1:-1]
                #else:
                #    raise ValueError('probably should return non endpoint intersection')
            elif (t0 == 1 and t1 > t0) or (t0 < t1 and t1 == 0):
                # if include single endpoint overlaps
                if endpoint:
                    return u.cp() if (u.isnear(p) or u.isnear(q)) else v.cp()
    elif uvcrspq == 0 and not upcrsuv == 0:
        return
    elif not uvcrspq == 0:
        if skew:
            upcrspq = near(utop.crs(ptoq).z, 0)
            t0 = upcrsuv / uvcrspq
            t1 = upcrspq / uvcrspq
            #t0 = near(near(t0, 0), 1)
            #t1 = near(near(t1, 0), 1)
            if u.isnear(p) or u.isnear(q) or v.isnear(p) or v.isnear(q):
                if endpoint or endpoints:
                    return u + utov * t1
            elif (0 <= t0 and t0 <= 1) and (0 <= t1 and t1 <= 1):
                return u + utov * t1


def bbox(points):
    a = points[0].cp()
    if len(points) == 1:
        b = a.cp()
    else:
        b = points[1].cp()
        for p in points:
            if p.x < a.x:
                a.x = p.x
            if p.y < a.y:
                a.y = p.y
            if p.z < a.z:
                a.z = p.z
            if p.x > b.x:
                b.x = p.x
            if p.y > b.y:
                b.y = p.y
            if p.z > b.z:
                b.z = p.z
    return a, b


def graham_scan(points):
    from .vec3 import vec3
    """compute convex hull in xy plane for a set of points"""
    # find the lowest y-coordinate and leftmost point, called P0
    i = 0
    for j in range(1, len(points)):
        if points[j].y <= points[i].y:
            #if (points[j].x < points[i].x) or (points[j].y < points[i].y):
            if (points[j].x <= points[i].x) or (points[j].y <= points[i].y):
                i = j
    p0 = points[i]
    # sort points by polar angle with P0,
    # ## if several points have the same polar angle then only keep the farthest
    x = vec3.X()
    angle = lambda j: 0 if i == j else x.saxy(points[j] - p0)
    order = sorted(range(len(points)), key=angle)
    # pop the last point from the stack if we turn clockwise to reach this point
    ccw = lambda u, v, w: (v - u).crs(w - v).z
    stack = []
    for j in order:
        while len(stack) > 1 and ccw(stack[-2], stack[-1], points[j]) < 0:
            stack.pop(-1)
        stack.append(points[j])
    return stack


def circumcircle(u, v, w):
    """find the center/radius of circumcircle of triangle uvw"""
    vu, wv, uw = (u - v), (v - w), (w - u)
    d = 2 * ((u - v).crs(v - w)).dot((u - v).crs(v - w))
    a = (v - w).dot(v - w) * (u - v).dot(u - w) / d
    b = (u - w).dot(u - w) * (v - u).dot(v - w) / d
    c = (u - v).dot(u - v) * (w - u).dot(w - v) / d
    o = u * a + v * b + w * c
    r = (u - o).mag()
    return o, r


def subdivide_triangles(old):
    new = []
    for u, v, w in old:
        uv, vw, wu = u.lerp(v, 0.5), v.lerp(w, 0.5), w.lerp(u, 0.5)
        new.append(( u, uv, wu))
        new.append(( v, vw, uv))
        new.append(( w, wu, vw))
        new.append((uv, vw, wu))
    return new


def slide(loop, n=1, m=0):
    queue = loop[:] + loop[:n]
    while len(queue) > n + m:
        yield queue[:n]
        queue.pop(0)


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


def loopO(loop):
    """Return a lexicographical loop origin index."""
    from .vec3 import vec3

    N = loop_normal(loop)
    X = vec3.Z().crs(N).nrm()
    Y = N.crs(X).nrm()

    loop = [vec3(p.dot(X), p.dot(Y), 0) for p in loop]

    sortkey = lambda l: (l[1].z, l[1].x, l[1].y)
    O = sorted(enumerate(loop), key=sortkey)[0][0]
    return O

    #zs = [p.z for p in loop]
    #return zs.index(min(zs))


def loop_area(loop):
    area = 0.0
    for i in range(len(loop)):
        u, v = loop[i - 1], loop[i]
        area -= (u.x + v.x) * (u.y - v.y) / 2.0
    return area


def loop_normal(loop):
    from .vec3 import vec3
    pn = vec3.O()
    for u, v, w in slide(loop, 3):
        uv, vw = (v - u), (w - v)
        #alpha = uv.ang(vw)
        #if alpha > 0:
        #    pn.trn(uv.crs(vw).nrm())
        pn.trn(uv.crs(vw))
    assert not pn.isO(), 'zero vector in loop_normal!'
    return pn.nrm()


def loop_contains(loop, other):
    #raise NotImplementedError('NOTE: current implementation is not reliable')
    #print('NOTE: current implementation is not reliable')
    for p in loop:
        if p.inbxy(other, True):
            return False
    for q in other:
        if q.inbxy(loop, True):
            return True
    else:
        return False


def loop_contains_triangle(loop, a, b, c):
    if not (a.inbxy(loop) or a.onbxy(loop, ie=True)):
        return False
    elif not (b.inbxy(loop) or b.onbxy(loop, ie=True)):
        return False
    elif not (c.inbxy(loop) or c.onbxy(loop, ie=True)):
        return False
    else:
        from .vec3 import vec3
        return vec3.com((a, b, c)).inbxy(loop)


def loop_exterior(loops):
    """find the one loop in loops which contains the others"""
    n_loops = len(loops)
    exterior = 0
    if n_loops > 1:
        for i in range(1, n_loops):
            if loop_contains(loops[i], loops[exterior]):
                exterior = i
    return exterior


def loop_offset(loop, r, closed=True, r0=10000):
    """loop contract for simple cases"""
    from .vec3 import vec3
    from .quat import quat

    n = loop_normal(loop)
    q = quat.toxy(n)
    q.rot(loop)

    Z = vec3.Z()
    offset = []
    for x, y, z in slide(loop, 3):
        #a = (y - x).ang(z - y)
        a = (y - x).axy(z - y)
        if isnear(a, 0, epsilon=0.01):
            p = y + ((z - x).crs(Z).nrm() * -r)
            offset.append(p)
        elif isnear(a, np.pi, epsilon=0.01):
            p = y + ((z - y).crs(Z).nrm() * r)
            offset.append(p)
            p = y + ((y - z).crs(Z).nrm() * r)
            offset.append(p)
        else:
            xy = (y - x).nrm()
            yz = (z - y).nrm()
            u = xy.crs(Z)
            v = yz.crs(Z)
            a = x + (u * -r) + (xy * -r0)
            b = y + (u * -r) + (xy *  r0)
            c = y + (v * -r) + (yz * -r0)
            d = z + (v * -r) + (yz *  r0)
            p = sintsxyp(a, b, c, d)
            offset.append(p)
    offset.insert(0, offset.pop(-1))

    if not closed:
        offset[ 0] = loop[ 0] + Z.crs(loop[ 1] - loop[ 0]).nrm() * r
        offset[-1] = loop[-1] + Z.crs(loop[-1] - loop[-2]).nrm() * r

    q = q.fp()
    q.rot(offset)
    q.rot(loop)

    return offset


def loop_contract(loop, r):
    from .vec3 import vec3
    from .mesh import planargraph
    if loop_normal(loop).z < 0:
        loop.reverse()
    newloop = []
    for i in range(len(loop)):
        u, v, w = loop[i - 2], loop[i - 1], loop[i]
        uv, vw = u.tov(v).nrm(), v.tov(w).nrm()
        uvn, vwn = uv.crs(vec3.Z()), vw.crs(vec3.Z())
        if u.isnear(w):
            newloop.append(v - uvn * r)
            newloop.append(v - vwn * r)
        elif uvn.isnear(vwn):
            newloop.append(v - uvn * r)
        else:
            a, b = u - uvn * r - uv * 1000, v - uvn * r + uv * 1000
            c, d = v - vwn * r - vw * 1000, w - vwn * r + vw * 1000
            ip = sintsxyp(a, b, c, d, False, False, False, True)
            assert ip is not None
            newloop.append(ip)
    segs = [(newloop[i - 1], newloop[i]) for i in range(len(newloop))]
    pg = planargraph(segs)
    pg.destem()
    loops = [[pg.vertices[i] for i in l] for l in pg.loops()]
    for loop in loops:
        if loop_normal(loop).z < 0:
            loop.reverse()
    loops = sorted(loops, key=lambda l: abs(loop_area(l)), reverse=True)
    return loops[1]


def loop_split(loop, maxlen):
    """Return new loop where no edge is greater than maxlen in length.

    Args:
        loop (seq): Sequence of vec3 instances.
        maxlen (float): Maximum edge length in resulting loop.

    Returns:
        list of new vec3 instances forming new loop.

    """
    newloop = []
    for i in range(len(loop)):
        u, v = loop[i - 1], loop[i]
        newloop.append(u.cp())
        el = u.tov(v).mag()
        if el > maxlen:
            n = 1
            while el / n > maxlen:
                n += 1
            newloop.extend(u.line(v, n - 1))
    return newloop
def edge_split(loop, maxlen=10):
    raise NotImplementedError('is this trash?')
    newloop = []
    for i in range(1, len(loop)):
        u, v = loop[i - 1], loop[i]
        newloop.append(u.cp())
        el = u.tov(v).mag()
        if el > maxlen:
            n = 1
            while el / n > maxlen:
                n += 1
            newloop.extend(u.line(v, n - 1))
    return newloop


def loop_smooth(loop, weight=0.1, iterations=1):
    from .vec3 import vec3
    for j in range(iterations):
        dps = []
        for i in range(len(loop)):
            u, v, w = loop[i - 2], loop[i - 1], loop[i]
            dps.append((vec3.com((u, w)) - v) * weight)
        dps.append(dps.pop(0))
        loop = [(p + dp) for p, dp in zip(loop, dps)]
    return loop


