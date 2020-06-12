"""View a 3D mesh of SRTM1 data in PNG format at various mesh resolutions"""
from meshmaker.plt import *
from meshmaker.mgl import show
from meshmaker.vec3 import vec3
from meshmaker.mesh import Mesh
from meshmaker.model import Model
from meshmaker.tform import TForm
from meshmaker.field import scalar_field
from pyglet.window import key
import cv2
import argparse


class ProgressiveTerrain:
    """Class which represents grayscale heightmaps with a mesh
    that can be progressively subdivided for better resolution"""

    @staticmethod
    def height_field(img):
        pixels = len(img)
        origin = vec3(pixels / 2, pixels / 2, 0)
        radius = pixels // 2
        deltaz = radius
        wtoi = TForm(t=vec3(pixels / 2, pixels / 2, 0) - origin,
                     s=vec3(pixels / (2 * radius), pixels / (2 * radius), 1))
        field = scalar_field.topography(img, wtoi, 1.0)
        field.origin = origin
        field.radius = radius
        field.deltaz = deltaz
        field.z0, field.dz = img.min(), img.max() - img.min()
        return field

    def __init__(self, height_img, sealevel=0.25):
        self.height = self.height_field(height_img)
        self.sealevel = sealevel
        self.divisions = 0

        self.dlevel = 0.1
        self.parameters = {
            key.RIGHT: self.divide,
            key.LEFT: self.undivide,
            key.UP: self.increase_sealevel,
            key.DOWN: self.decrease_sealevel,
            key.K: self.increase_deltaz,
            key.J: self.decrease_deltaz,
        }

        print('Controls:')
        for button, func in self.parameters.items():
            name = func.__name__[func.__name__.rfind('.'):]
            print(f'Key: {button} -> {name}')

    def divide(self, *ags):
        self.divisions += 1

    def undivide(self, *ags):
        self.divisions = max(0, self.divisions - 1)

    def increase_deltaz(self, *ags):
        self.height.deltaz *= 1.1

    def decrease_deltaz(self, *ags):
        self.height.deltaz /= 1.1

    def increase_sealevel(self, *ags):
        self.sealevel += 0.05

    def decrease_sealevel(self, *ags):
        self.sealevel -= 0.05

    def scene(self):
        """Remake the land/sea meshes and return them in a scenegraph
        (NOTE: This method is called by the rendering window)"""
        height = self.height
        land, sea = Mesh(), Mesh()
        for u, v, w in height.origin.fan(height.radius, 6, True):
            sea.af((u.cp(), v.cp(), w.cp()))
            land.af((u.cp(), v.cp(), w.cp()))
        for i in range(self.divisions):
            land = land.subdivide()
        print(f'Triangles: {len(land.faces) + len(sea.faces)}')
        sealevel = (height.z0 + self.sealevel * height.dz)
        for v in sea.vertices:
            if v is not None:
                v.z = sealevel * height.deltaz
        for v in land.vertices:
            if v is not None:
                v.z = height(v) * height.deltaz
        land.uvs = land.project_uvs_xy(0.1)
        sea.uvs = sea.project_uvs_xy(0.1)
        land.normals = land.vertex_normals(smooth=True)
        root = TForm(s=vec3.U(20.0 / 3600.0))
        root.add(TForm(models=[Model(meshes={'generic_13': [land]})]))
        root.add(TForm(models=[Model(meshes={'generic_10': [sea]})]))
        return root


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='View SRTM1 data in PNG format')
    parser.add_argument('--png', type=str,
                        default='./srtm1/N38W112.png',
                        help='Which PNG to visualize')
    args = parser.parse_args()

    height = cv2.imread(args.png, cv2.IMREAD_GRAYSCALE)
    height = height / height.max()
    height = height[400:600, 400:600]

    f, ax = plot()
    im = ax.imshow(height, origin='lower')
    plt.show()

    show(ProgressiveTerrain(height),
         texture_directory='../../resources/textures')
