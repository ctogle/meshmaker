import matplotlib.pyplot as plt


def plot(i=1, j=1, figsize=(8, 8)):
    f, ax = plt.subplots(i, j, figsize=figsize, dpi=100)
    return f, ax


def annotate_point(ax, p, annotation):
    ax.annotate(annotation, xy=(p.x, p.y), xytext=(-10, 10),
        textcoords='offset points', ha='right', va='bottom',
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    return ax


def plot_point(ax, p, col='k', mk='o', annotation=None):
    ax.plot(p.x, p.y, color=col, marker=mk)
    if annotation is not None:
        ax = annotate_point(ax, p, annotation)
    return ax


def plot_edge(ax, u, v, lw=3, ls='-', col='k'):
    ax.plot([u.x, v.x], [u.y, v.y], lw=lw, ls=ls, color=col)
    return ax


def plot_edges(ax, edges, lw=3, ls='-', col='k', mk=None):
    edges = [(edges[i - 1], edges[i]) for i in range(1, len(edges))]
    for i, (u, v) in enumerate(edges):
        if mk:
            #plot_point(ax, u, col, mk=mk, annotation=f'{i}')
            plot_point(ax, v, col, mk=mk, annotation=f'{i}')
        plot_edge(ax, u, v, lw, ls, col)
    return ax


def plot_loop(ax, loop, lw=3, ls='-', col='k', mk=None):
    return plot_edges(ax, [loop[-1]] + loop, lw, ls, col, mk)


def plot_pg(ax, pg, lw=3, ls='-', col='k', mk='o', annotate=True):
    for i, v in enumerate(pg.vertices):
        if v is None:
            continue
        ax.plot([v.x], [v.y], marker=mk, color=col)
        if annotate:
            text = v.properties.get('annotation')
            if text is None:
                text = str(pg.rings[i]) if len(pg.rings[i]) > 2 else ''
                #text = str(pg.rings[i]) if len(pg.rings[i]) > -1 else ''
            if text:
                text = '%d\n%s' % (i, text)
                annotate_point(ax, v, text)
    for i, j, properties in (e for e in pg.edges if e is not None):
        u, v = pg.vertices[i], pg.vertices[j]
        if u is not None and v is not None:
            plot_edge(ax, u, v, lw=lw, ls=ls, col=col)
    return ax
