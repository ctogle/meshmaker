# subtract one mesh from another
from meshmaker.mesh import Mesh
from meshmaker.vec3 import vec3
from meshmaker.mgl import show

a = vec3.X(-1).ring(3, 4)
a = Mesh.prism(a, 2, 1.0)

b = vec3.X(+1).ring(2, 4)
b = Mesh.prism(b, 2, 2.0)

c = a.difference(b)
show(c)
