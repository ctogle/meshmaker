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
            text = '%d\n%s' % (i, str(pg.rings[i]))
            annotate_point(ax, v, text)
    for i, j in filter(lambda e: bool(e), pg.edges):
        u, v = pg.vertices[i], pg.vertices[j]
        ax.plot([u.x, v.x], [u.y, v.y], ls=ls, lw=lw, color=col)
    return ax


def plot_basis(ax, b, p, which=None):
    l = 1.5
    major, minor = b.evectors()
    if (which == 'major' or not which) and major:
        dp = major.nrm() * l
        plot_point(ax, p + dp, col='b')
        plot_edge(ax, p, p + dp, col='b', lw=2)
    if (which == 'minor' or not which) and minor:
        dp = minor.nrm() * l
        plot_point(ax, p + dp, col='g')
        plot_edge(ax, p, p + dp, col='g', lw=2)
    plot_point(ax, p, col='r')
    return ax


def plot_field(ax, f, ps, which=None):
    for p in ps:
        t = f.basis(p)
        plot_basis(ax, t, p, which=which)
    return ax


