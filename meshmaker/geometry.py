

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


def orient2d(a, b, c):
    """Signed area of the triangle a - c, b - c
        0 if a, b, and c are colinear
    """
    m11 = a.x - c.x
    m12 = a.y - c.y
    m21 = b.x - c.x
    m22 = b.y - c.y
    det = m11 * m22 - m12 * m21
    return near(det, 0)


def sintsxyp(u, v, p, q, endpoint=True, endpoints=True, colinear=True, skew=True):
    utov, ptoq, utop = u.tov(v), p.tov(q), u.tov(p)
    uvcrspq, upcrsuv = near(utov.crs(ptoq).z, 0), near(utop.crs(utov).z, 0)
    assert(isnear(utov.z, 0) and isnear(ptoq.z, 0))
    if uvcrspq == 0 and upcrsuv == 0:
        if colinear:
            uv, pq = utov.mag(), ptoq.mag()
            uvdotpq = uv.dot(pq), up.dot(uv)
            updotuv = up.dot(uv)
            assert(not (uvdotpq == 0))
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
            if ((t0 == 0 or t0 == 1) and (t1 == 0 or t1 == 1)):
                if endpoints:
                    return u + utov * t1
            elif not endpoint and ((t0 == 0 or t1 == 0) or (t0 == 1 or t1 == 1)):
                if endpoint:
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
            if p.x > b.x:
                b.x = p.x
            if p.y > b.y:
                b.y = p.y
    return a, b


def subdivide_triangles(old):
    new = []
    for u, v, w in old:
        uv, vw, wu = u.lerp(v, 0.5), v.lerp(w, 0.5), w.lerp(u, 0.5)
        new.append(( u, uv, wu))
        new.append(( v, vw, uv))
        new.append(( w, wu, vw))
        new.append((uv, vw, wu))
    return new


def slide(loop, n=1):
    queue = loop[:] + loop[:n]
    while len(queue) > n:
        yield queue[:n]
        queue.pop(0)


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
        alpha = uv.axy(vw)
        pn.trn(uv.crs(vw).nrm())
    return pn.nrm()


def loop_contains(loop, other):
    for p in loop:
        if p.inbxy(other, True):
            return False
    for q in other:
        if q.inbxy(loop, True):
            return True
    else:
        return False


def loop_exterior(loops):
    """find the one loop in loops which contains the others"""
    n_loops = len(loops)
    exterior = 0
    if n_loops > 1:
        for i in range(1, n_loops):
            if loop_contains(loops[i], loops[exterior]):
                exterior = i
    return exterior


def loop_contract(loop, r):
    from .vec3 import vec3
    from .mesh import planargraph
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
    loops = sorted(loops, key=lambda l: abs(loop_area(l)), reverse=True)
    return loops[1]


def loop_split(loop, maxlen=10):
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


def loop_smooth(loop, weight=0.1):
    from .vec3 import vec3
    dps = []
    for i in range(len(loop)):
        u, v, w = loop[i - 2], loop[i - 1], loop[i]
        dps.append((vec3.com((u, w)) - v) * weight)
    dps.append(dps.pop(0))
    return [(p + dp) for p, dp in zip(loop, dps)]


if __name__ == '__main__':
    from .vec3 import vec3

    (u, v, q, p) =\
        (vec3(0.0000, 2.0000, 0.0000),
         vec3(0.4337, 2.9011, 0.0000),
         vec3(-0.4337, 1.9011, 0.0000),
         vec3(0.0000, 1.0000, 0.0000))
    print(u, v, q, p)

    ip = sintsxyp(u, v, q, p, False, False, False, True)
    print(ip)

    from .plt import plot, plot_pg, plot_point, plot_edge
    f, ax = plot(figsize=(10, 10))
    plot_edge(ax, u, v, col='g', lw=5)
    plot_edge(ax, p, q, col='b', lw=5)
    plot_point(ax, u, col='g', mk='s', annotation='u')
    plot_point(ax, v, col='g', mk='s', annotation='v')
    plot_point(ax, p, col='b', mk='s', annotation='p')
    plot_point(ax, q, col='b', mk='s', annotation='q')
    plot_point(ax, ip, col='r', mk='s', annotation='ip')
    import matplotlib.pyplot as plt
    plt.show()

