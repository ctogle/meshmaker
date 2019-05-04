from io import StringIO
import os


mtl = '''
newmtl {0}
Ka 1.0 1.0 1.0
Kd 1.0 1.0 1.0
Ks 0.0 0.0 0.0
d 1.0
illum 2
map_Kd {1}
'''

def mtl_file(textures):
    f = StringIO()
    f.write('# Material Count: {0}\n'.format(len(textures)))
    for texture in textures:
        f.write(mtl.format(texture, textures[texture]))
    return f.getvalue()


def obj_file(model, mtl_path):
    f = StringIO()
    f.write('# your header here\n')
    f.write('mtllib %s\n' % mtl_path)
    f.write('o %s\n' % model.name)
    for mtl in model.surfaces:
        for surface in model.surfaces[mtl]:
            for p, n, u in surface.vertices:
                f.write( 'v %0.6f %0.6f %0.6f\n' % (p.x, p.y, p.z))
                f.write('vn %0.6f %0.6f %0.6f\n' % (n.x, n.y, n.z))
                f.write('vt %0.6f %0.6f\n'       % (u.x, u.y))
            f.write('usemtl %s\n' % mtl)
            f.write('s off\n')
            for v1, v2, v3 in surface.faces:
                face = (v1 + 1, v2 + 1, v3 + 1)
                f.write('f {0}/{0}/{0} {1}/{1}/{1} {2}/{2}/{2}\n'.format(*face))
    return f.getvalue()


def obj_world(prefix, models, textures):
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    mtl_path = '{0}materials.mtl'.format(prefix)
    with open(mtl_path, 'w') as f:
        f.write(mtl_file(textures))
    for m in models:
        obj_path = '{0}{1}.obj'.format(prefix, m.name)
        with open(obj_path, 'w') as f:
            f.write(obj_file(m, mtl_path))


