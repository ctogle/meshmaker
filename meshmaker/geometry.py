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


def slide(iterable, n=1, m=0):
    queue = iterable[:] + iterable[:n]
    while len(queue) > n + m:
        yield queue[:n]
        queue.pop(0)


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]






# TODO: remove
def loop_area(loop):
    area = 0.0
    for i in range(len(loop)):
        u, v = loop[i - 1], loop[i]
        area -= (u.x + v.x) * (u.y - v.y) / 2.0
    return area


# TODO: remove
def loop_normal(loop):
    from .vec3 import vec3
    pn = vec3.O()
    for u, v, w in slide(loop, 3):
        uv, vw = (v - u), (w - v)
        pn.trn(uv.crs(vw))
    assert not pn.isO(), 'zero vector in loop_normal!'
    return pn.nrm()
