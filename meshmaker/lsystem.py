from .vec3 import vec3
from .quat import quat
from .mesh import planargraph
from io import StringIO

import numpy as np
import matplotlib.pyplot as plt


def lsystem(axiom, rules, position, direction, iterations, dtheta):
    state = axiom
    for i in range(iterations):
        out = StringIO()
        for c in state:
            if c in rules:
                out.write(rules[c])
            else:out.write(c)
        state = out.getvalue()
    #pg = planargraph()
    segs = []
    tip = (position, direction, None)
    for c in state:
        if c == '{':
            tip = (tip[0].cp(), tip[1].cp(), tip)
        elif c == '}':
            tip = tip[2]
        elif c == '-':
            #pg.ae(tip[0].cp(), tip[0].cp() + tip[1].nrm(), 0.1)
            segs.append((tip[0].cp(), tip[0].cp() + tip[1].nrm()))
            tip[0].trn(tip[1].nrm())
        elif c == '[':
            tip[1].rot(quat.av( dtheta, vec3(0, 0, 1)))
        elif c == ']':
            tip[1].rot(quat.av(-dtheta, vec3(0, 0, 1)))
    #return pg
    return segs


if __name__ == '__main__':
    dragon_curve = (
        '-X',
        {'X': 'X[Y-[', 'Y': ']-X]Y'},
        vec3(0, 0, 0),
        vec3(0, 1, 0),
        5,
        np.pi / 2.0,
    )
    #pg = lsystem(*dragon_curve)

    pythagoras_tree = (
        'X',
        {'-': '--', 'X': '-{[X}]X'},
        vec3(0, 0, 0),
        vec3(1, 0, 0),
        5,
        np.pi / 12.0,
    )
    #pg = lsystem(*pythagoras_tree)

    axial_tree = (
       'X',
        {'X': '-{[X}{]X}-X', '-': '--'},
        vec3(0, 0, 0),
        vec3(0, 1, 0),
        5,
        25.7 * (np.pi / 180),
    )
    #pg = lsystem(*axial_tree)

    plant = (
       'X',
        {'X': '-]{{X}[X}[-{[-X}]X', '-': '--'},
        vec3(0, 0, 0),
        vec3(0, 1, 0),
        3,
        25 * (np.pi / 180),
    )
    #pg = lsystem(*plant)

    grass = (
       '-',
        {'-': '-{[-}-{]-}-'},
        vec3(0, 0, 0),
        vec3(0, 1, 0),
        2,
        25.7 * (np.pi / 180),
    )
    pg = planargraph(lsystem(*grass))

    from .plt import plot, plot_pg
    f, ax = plot(figsize=(10, 10))
    plot_pg(ax, pg, annotate=False)
    ax.set_title('%d vertices, %d edges' % (pg.vertex_count, pg.edge_count))
    plt.show()


