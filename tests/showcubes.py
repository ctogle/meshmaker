from meshmaker.model import Model
from meshmaker.mesh import Mesh
from meshmaker.tform import TForm
from meshmaker.vec3 import vec3
from meshmaker.quat import quat
from meshmaker.mgl import show, MainShader, EdgeShader
from functools import partial
import numpy as np

# create a simpler calling method to show things via OpenGL
shaders = [MainShader(), EdgeShader()]
show = partial(show, programs=shaders, background=vec3.U(0.5))

# make a cube and some copies
a = Mesh.cube_mesh(r=1)
b = a.cp()
c = a.cp()
d = a.cp()

# pick some t, q, s (translation, rotation, scale)
t, q, s = vec3.X(2), quat.av(np.pi / 6, vec3.Z()), vec3(0.6, 0.8, 1.2)

# transform some cubes with t, q, s
t.trnps(
    q.rot(
        s.sclps(
            b.vertices)))
t.trnps(t.trnps(
    q.rot(q.rot(
        s.sclps(s.sclps(
            c.vertices))))))
t.trnps(t.trnps(t.trnps(
    q.rot(q.rot(q.rot(
        s.sclps(s.sclps(s.sclps(
            d.vertices)))))))))

# put cubes in a model in a scenegraph to assign a texture to each cube
meshes={
    'generic_8': [a],
    'generic_9': [b],
    'generic_10': [c],
    'generic_11': [d],
        }
show(TForm(t=vec3.X(-2), models=[Model(meshes=meshes)]))
