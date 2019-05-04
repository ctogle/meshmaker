from .mesh import trimesh, planargraph
from .model import model
from .vec3 import vec3
from .obj import obj_world

import matplotlib.pyplot as plt


def cube():
    p1 = vec3(0, 0, 0)
    p2 = vec3(1, 0, 0)
    p3 = vec3(1, 1, 0)
    p4 = vec3(0, 1, 0)
    p5 = vec3(0, 1, 1)
    p6 = vec3(1, 1, 1)
    p7 = vec3(1, 0, 1)
    p8 = vec3(0, 0, 1)
    m = trimesh()
    m.af(p4, p3, p2, p1)
    m.af(p8, p7, p6, p5)
    m.af(p1, p2, p7, p8)
    m.af(p2, p3, p6, p7)
    m.af(p3, p4, p5, p6)
    m.af(p4, p1, p8, p5)
    return m


def river():
    p1 = vec3(0, 0, 0)
    p2 = vec3(1, 0, 0)
    p3 = vec3(1, 1, 0)
    p4 = vec3(0, 1, 0)
    pg = planargraph()
    pg.ae(p1, p2)
    pg.ae(p2, p3)
    pg.ae(p3, p4)
    pg.ae(p4, p1)

    i, j = pg.edges[0]
    z = vec3(0, 0, 1)
    pg.follow(i, j, z)

    return pg


if __name__ == '__main__':
    r = river()

    f, ax = plt.subplots(1, 1, figsize=(8, 8))
    for v in r.vertices:
        ax.plot([v.x], [v.y], [v.z], marker='o', color='g')
    for i, j in r.edges:
        u, v = r.vertices[i], r.vertices[j]
        ax.plot([u.x, v.x], [u.y, v.y], [u.z, v.z],
                linestyle='-', linewidth=3, color='g')
    plt.show()

    quit()

    prefix = './obj.world/__'

    textures = {
        'default': '../resources/textures/orangeboxtex.png',
            }

    models = (
        model().add(cube(), 'default'),
        model().add(cube(), 'default'),
    )

    obj_world(prefix, models, textures)

