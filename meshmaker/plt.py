import matplotlib.pyplot as plt


def plot(i=1, j=1, figsize=(8, 8)):
    f, ax = plt.subplots(i, j, figsize=figsize)
    return f, ax


def annotate_point(ax, p, annotation):
    ax.annotate(annotation, xy=(p.x, p.y), xytext=(-10, 10),
        textcoords='offset points', ha='right', va='bottom',
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    return ax


def plot_point(ax, p, col='k', mk='o', annotation=None):
    ax.plot(p.x, p.y, color=col, marker=mk)
    if annotation:
        ax = annotate_point(ax, p, annotation)
    return ax


def plot_edge(ax, u, v, lw=3, ls='-', col='k'):
    ax.plot([u.x, v.x], [u.y, v.y], lw=lw, ls=ls, color=col)
    return ax


def plot_loop(ax, loop, lw=3, ls='-', col='k', mk=None):
    for i in range(len(loop)):
        if mk:
            plot_point(ax, loop[i], col, mk=mk, annotation=f'{i}')
        plot_edge(ax, loop[i - 1], loop[i], lw, ls, col)
    return ax


def plot_pg(ax, pg, lw=3, ls='-', col='k', mk='o', annotate=True):
    for i, v in enumerate(pg.vertices):
        if v is None:
            continue
        ax.plot([v.x], [v.y], marker=mk, color=col)
        if annotate:
            text = '%d\n%s' % (i, str(pg.rings[i]))
            annotate_point(ax, v, text)
    for i, j in filter(lambda e: bool(e), pg.edges):
        u, v = pg.vertices[i], pg.vertices[j]
        ax.plot([u.x, v.x], [u.y, v.y], ls=ls, lw=lw, color=col)
    return ax


