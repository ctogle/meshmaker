import os
import numpy as np
from pyrr import Matrix44
from PIL import Image
import moderngl
import moderngl_window
from moderngl_window.conf import settings
from moderngl_window.timers.clock import Timer
from .vec3 import vec3
from .base import Base, Laziness
from .tform import TForm
from .model import Model
from .geometry import bbox


class LazyMaterials(Laziness):

    def __init__(self, ctx, texture_directory, **kws):
        self.textures = {}
        for r, ds, fs in os.walk(texture_directory):
            for f in fs:
                name = f[:f.rfind('.')]
                path = os.path.join(r, f)
                self.textures[name] = path
        super().__init__(self.method(ctx), **kws)

    def __iter__(self):
        yield from self.textures.keys()

    def method(self, ctx):
        def wrapped(name):
            path = self.textures[name]
            img = Image.open(path)
            imgbytes = img.tobytes()
            bytesize = int(len(imgbytes) / np.product(img.size))
            texture = ctx.texture(img.size, bytesize, imgbytes)
            print(f'Loaded texture: {name} ({path})')
            return texture
        return wrapped


class Camera(Base):

    def pan(self, dx, dy):
        look, left, up = self.basis
        self.target += (left * dx + up * dy) * self.speed

    def mouse_drag_event(self, x, y, dx, dy):
        if self._last_mouse == 1:
            self.pan(dx, dy)
        elif self._last_mouse == 2:
            self.theta -= dx * self.speed
            self.phi -= dy * self.speed

    def mouse_scroll_event(self, x_offset, y_offset):
        self.radius -= y_offset * self.speed * 10

    def mouse_press_event(self, x, y, button):
        self._last_mouse = button

    def key_event(self, key, action, modifiers):
        if action == self.keys.ACTION_PRESS:
            if key == self.keys.R:
                self.reset_view()
            elif key == self.keys.W:
                self.pan(  0, 10)
            elif key == self.keys.S:
                self.pan(  0,-10)
            elif key == self.keys.A:
                self.pan( 10,  0)
            elif key == self.keys.D:
                self.pan(-10,  0)

    def reset_view(self):
        self.target = self.default_target.cp()
        self.radius = self.default_radius
        self.theta  = -np.pi / 3
        self.phi    = np.pi / 3
        self.aspect_ratio = 16 / 9
        self.speed = 0.01

    def __init__(self, keys, prog, radius=None, target=None):
        self.default_radius = 10 if radius is None else radius
        self.default_target = vec3.O() if target is None else target
        self.keys = keys
        self.mvp = prog['Mvp']
        self.light = prog['Light']
        self.reset_view()

    @property
    def xyz(self):
        x = self.radius * np.cos(self.theta) * np.sin(self.phi)
        y = self.radius * np.sin(self.theta) * np.sin(self.phi)
        z = self.radius * np.cos(self.phi)
        return x + self.target.x, y + self.target.y, z + self.target.z

    @property
    def basis(self):
        look = (self.target - vec3(*self.xyz)).nrm()
        left = vec3.Z().crs(look)
        up = look.crs(left)
        return look, left, up

    def __call__(self):
        x, y, z = self.xyz
        proj = Matrix44.perspective_projection(
            45.0, self.aspect_ratio, 0.1, 100.0)
        lookat = Matrix44.look_at(
            (            x,             y,             z),
            (self.target.x, self.target.y, self.target.z),
            (            0,             0,             1))
        self.mvp.write((proj * lookat).astype('f4').tobytes())
        self.light.value = (x, y, z)


class Window(Base):

    shaders = {
        'vertex_shader': '''
            #version 330

            uniform mat4 Mvp;

            in vec3 in_position;
            in vec3 in_normal;
            in vec2 in_texcoord_0;

            out vec3 v_vert;
            out vec3 v_norm;
            out vec2 v_text;

            void main() {
                gl_Position = Mvp * vec4(in_position, 1.0);
                v_vert = in_position;
                v_norm = in_normal;
                v_text = in_texcoord_0;
            }
        ''',
        'fragment_shader': '''
            #version 330

            uniform vec3 Light;
            uniform sampler2D Texture;

            in vec3 v_vert;
            in vec3 v_norm;
            in vec2 v_text;

            out vec4 f_color;

            void main() {
                float lum = clamp(dot(normalize(Light - v_vert), normalize(v_norm)), 0.0, 1.0) * 0.3 + 0.4;
                vec3 base = vec3(0.5, 0.5, 0.5) * lum;
                vec3 spec = vec3(1.0, 1.0, 1.0) * pow(lum, 5.7);
                vec4 tex = texture(Texture, v_text);
                f_color = vec4(base * 0.1 + tex.rgb * lum + spec, tex.a);
            }
        ''',
    }

    def __init__(self, source, texture_directory, **kws):
        self.source = source
        self.targets = []

        settings.WINDOW['class'] = 'moderngl_window.context.pyglet.Window'
        self.wnd = moderngl_window.create_window_from_settings()

        self.ctx = self.wnd.ctx
        self.ctx.enable(moderngl.DEPTH_TEST | moderngl.CULL_FACE)

        self.prog = self.ctx.program(**self.shaders)

        self.materials = LazyMaterials(self.ctx, texture_directory)

        self.camera = Camera(self.wnd.keys, self.prog, **kws)

        self.wnd.mouse_drag_event_func = self.camera.mouse_drag_event
        self.wnd.mouse_scroll_event_func = self.camera.mouse_scroll_event
        self.wnd.mouse_press_event_func = self.camera.mouse_press_event
        self.wnd.key_event_func = self.key_event

        self.run()

    def key_event(self, key, action, modifiers):
        self.camera.key_event(key, action, modifiers)
        keys = self.wnd.keys
        if action == keys.ACTION_PRESS:
            if key == keys.Q:
                self.wnd.close()
            elif key == keys.U:
                self.update()
            elif key == keys.TAB:
                self.ctx.wireframe = not self.ctx.wireframe
            else:
                parameters = getattr(self.source, 'parameters', {})
                operation = parameters.get(key)
                if operation is not None:
                    operation(key, action, modifiers)
                    self.update()

    def run(self):
        timer = Timer()
        timer.start()
        self.update()
        while not self.wnd.is_closing:
            self.wnd.clear()
            time, frame_time = timer.next_frame()
            self.render(time, frame_time)
            self.wnd.swap_buffers()
        self.wnd.destroy()

    def render(self, time, frame_time):
        self.ctx.clear(0.9, 0.9, 1.0)
        self.ctx.enable(moderngl.DEPTH_TEST)
        self.camera()
        for mode, texture, vao in self.targets:
            if texture is not None:
                texture.use()
            vao.render(mode)

    def axes(self):
        # TODO: draw on top and with colors...
        axes = np.array([[[0, 0, 0], [1, 0, 0], [0, 0, 0]],
                         [[0, 0, 0], [0, 1, 0], [0, 0, 0]],
                         [[0, 0, 0], [0, 0, 1], [0, 0, 0]]])
        vbo = self.ctx.buffer(axes.astype('f4').tobytes())
        vao = self.ctx.simple_vertex_array(self.prog, vbo, 'in_position')
        return (moderngl.LINES, None, vao)

    def draw_material(self, tf, material, vertices=None, wireframe=None):
        """Recursively aggregate data from a scene graph for a material"""
        vertices = [] if vertices is None else vertices
        wireframe = [] if wireframe is None else wireframe
        if hasattr(tf, 'models'):
            for model in tf.models:
                for mesh in model.meshes.get(material, ()):
                    try:
                        tridata, wiredata = mesh.T2F_N3F_V3F(tf)
                        vertices.extend(tridata)
                        wireframe.extend(wiredata)
                    except:
                        print(f'failed to generate mesh {mesh}')
                        raise
        for child in tf.children:
            self.draw_material(child, material, vertices, wireframe)
        return vertices, wireframe

    def update(self):
        """Recompute the scenegraph and vao data (i.e. self.targets)"""
        if hasattr(self.source, 'scene'):
            self.scene = self.source.scene()
        else:
            self.scene = TForm()
            print(f'Source has no "scene" method: {self.source}')
        self.targets = [self.axes()]
        for name in self.materials:
            vertices, wireframe = self.draw_material(self.scene, name)
            if vertices:
                vertices = np.array(vertices)
                vbo = self.ctx.buffer(vertices.astype('f4').tobytes())
                vao = self.ctx.simple_vertex_array(self.prog, vbo,
                    'in_texcoord_0', 'in_normal', 'in_position')
                self.targets.append((moderngl.TRIANGLES, self.materials[name], vao))

    @classmethod
    def test(cls, texture_directory, **kws):
        cube = Model.cube_model('generic_0', 1)
        source = Base(scene=(lambda : TForm(models=[cube])))
        cls(source, texture_directory, **kws)
