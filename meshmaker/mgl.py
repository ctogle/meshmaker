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
from .mesh import Mesh
from .geometry import bbox, batch, slide


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

    def __init__(self, keys, radius=None, target=None, **kws):
        self.default_radius = 10 if radius is None else radius
        self.default_target = vec3.O() if target is None else target
        self.keys = keys
        self.reset_view()

        self.projection = Matrix44.perspective_projection(
                        45.0, self.aspect_ratio, 0.1, 100.0)
        #self.projection = Matrix44.orthogonal_projection(
        #                0.0, 1024, 768, 0, -1.0, 1.0)

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

    @property
    def view(self):
        x, y, z = self.xyz
        lookat = Matrix44.look_at(
            (            x,             y,             z),
            (self.target.x, self.target.y, self.target.z),
            (            0,             0,             1))
        return lookat

    @property
    def viewdir(self):
        return (self.target - vec3(*self.xyz)).nrm()


class ShaderProgram(Base):

    def _readiffile(self, s):
        if s is not None:
            return open(s, 'r').read().strip() if os.path.isfile(s) else s

    def __init__(self, vs=None, gs=None, fs=None,
                 signature=None, primitive=None, active=True):
        self.vs = vs
        self.gs = gs
        self.fs = fs
        self.signature = signature
        self.primitive = primitive
        self.active = active

    def build(self, ctx):
        """Compute shader program from context and associated GLSL scripts"""
        vs = self._readiffile(self.vs)
        gs = self._readiffile(self.gs)
        fs = self._readiffile(self.fs)
        self._ctx = ctx
        self._program = ctx.program(
            vertex_shader=vs,
            geometry_shader=gs,
            fragment_shader=fs)
        return self

    def init(self):
        """Set GLSL uniform variables at start"""
        pass

    def update(self, camera):
        """Set GLSL uniform variables per frame for camera"""
        pass

    def render(self, vertexdata, wiredata, material):
        """Compute function which renders vertexdata using shaders"""

        def f():
            pass

        return f


class MainShader(ShaderProgram):

    def __init__(self):
        vs = os.path.join(os.path.dirname(__file__), 'shaders', 'main.vs')
        fs = os.path.join(os.path.dirname(__file__), 'shaders', 'main.fs')
        signature = ('in_texcoord_0', 'in_normal', 'in_position')
        primitive = moderngl.TRIANGLES
        super().__init__(vs=vs, fs=fs, signature=signature, primitive=primitive)
        self.ka, self.kd, self.ks = 0.1, 0.5, 0.0

    def init(self):
        """Set GLSL uniform variables at start"""
        self._program["dirLight1.direction"].value = (-2.0,-1.0,-3.0)
        self._program["dirLight1.diffuse"].value   = (self.kd, self.kd, self.kd)
        self._program["dirLight1.ambient"].value   = (self.ka, self.ka, self.ka)
        self._program["dirLight1.specular"].value  = (self.ks, self.ks, self.ks)
        self._program["dirLight2.direction"].value = ( 2.0, 1.0, 3.0)
        self._program["dirLight2.diffuse"].value   = (self.kd, self.kd, self.kd)
        self._program["dirLight2.ambient"].value   = (self.ka, self.ka, self.ka)
        self._program["dirLight2.specular"].value  = (self.ks, self.ks, self.ks)
        self._program["material.specular"].value  = ( 0.9, 0.9, 0.9)
        self._program["material.shininess"].value = 100

    def update(self, camera):
        """Set GLSL uniform variables per frame for camera"""
        self._program['projection'].write(camera.projection.astype('f4').tobytes())
        self._program['view'].write(camera.view.astype('f4').tobytes())
        self._program['viewDir'].value = tuple(camera.viewdir)

    def render(self, vertexdata, wiredata, material):
        """Compute function which renders vertexdata using shaders"""
        if vertexdata:
            vertexdata = np.array(vertexdata)
            vbo = self._ctx.buffer(vertexdata.astype('f4').tobytes())
            vao = self._ctx.simple_vertex_array(self._program, vbo, *self.signature)

            def f():
                if material:
                    material.use()
                vao.render(self.primitive)
        else:
            f = lambda: None

        return f


class EdgeShader(ShaderProgram):

    def __init__(self):
        vs='''
            #version 330
            uniform mat4 projection;
            uniform mat4 view;

            in vec3 in_vert;

            void main() {
                mat4 Mvp = projection * view;
                gl_Position = Mvp * vec4(in_vert, 1.0);
            }
        '''
        fs='''
            #version 330
            out vec4 f_color;
            void main() {
                f_color = vec4(0.1, 1.0, 0.1, 1.0);
            }
        '''
        super().__init__(vs=vs, fs=fs,
                         primitive=moderngl.LINES,
                         signature=('in_vert', ))

    def init(self):
        pass

    def update(self, camera):
        self._program['projection'].write(camera.projection.astype('f4').tobytes())
        self._program['view'].write(camera.view.astype('f4').tobytes())

    def wiredata(self, vertexdata, dn=0.002):
        positions = []
        for i, (v1, v2, v3) in enumerate(batch(list(batch(vertexdata, 8)), 3)):
            n1, n2, n3 = vec3(*v1[2:5]), vec3(*v2[2:5]), vec3(*v3[2:5])
            p1, p2, p3 = vec3(*v1[5:8]), vec3(*v2[5:8]), vec3(*v3[5:8])
            p1.trn(n1 * dn)
            p2.trn(n2 * dn)
            p3.trn(n3 * dn)
            positions.extend(p1)
            positions.extend(p2)
            positions.extend(p2)
            positions.extend(p3)
            positions.extend(p3)
            positions.extend(p1)
        positions = np.array(positions)
        return positions

    def render(self, vertexdata, wiredata, material):
        if vertexdata:
            positions = self.wiredata(vertexdata)
            vbo = self._ctx.buffer(positions.astype('f4').tobytes())
            vao = self._ctx.simple_vertex_array(self._program, vbo, *self.signature)

            def f():
                vao.render(self.primitive)
        else:
            f = lambda: None

        return f


class WireShader(ShaderProgram):

    def __init__(self, color=vec3.X()):
        vs='''
            #version 330
            uniform mat4 projection;
            uniform mat4 view;

            in vec3 in_vert;

            void main() {
                mat4 Mvp = projection * view;
                gl_Position = Mvp * vec4(in_vert, 1.0);
            }
        '''
        fs='''
            #version 330
            out vec4 f_color;
            void main() {
                f_color = vec4(%0.2f, %0.2f, %0.2f, 1.0);
            }
        ''' % tuple(color)
        super().__init__(vs=vs, fs=fs,
                         primitive=moderngl.LINES,
                         signature=('in_vert', ))

    def init(self):
        pass

    def update(self, camera):
        self._program['projection'].write(camera.projection.astype('f4').tobytes())
        self._program['view'].write(camera.view.astype('f4').tobytes())

    def render(self, vertexdata, wiredata, material):
        positions = []
        #for u, v in slide(wiredata, 2):
        for u, v in wiredata:
            positions.extend(u)
            positions.extend(v)

        if positions:
            positions = np.array(positions)
            vbo = self._ctx.buffer(positions.astype('f4').tobytes())
            vao = self._ctx.simple_vertex_array(self._program, vbo, *self.signature)

            def f():
                vao.render(self.primitive)
        else:
            f = lambda : None

        return f


class Window(Base):

    def toggleshader(self, shader):
        def toggle():
            shader.active = not shader.active
            self.update()
        return toggle

    def __init__(self, source, texture_directory, **kws):
        super().__init__(**kws)
        self.source = source
        self.targets = []

        if not hasattr(self, 'background'):
            self.background = vec3(0.1, 0.1, 0.1)

        settings.WINDOW['class'] = 'moderngl_window.context.pyglet.Window'
        self.wnd = moderngl_window.create_window_from_settings()

        self.ctx = self.wnd.ctx
        self.ctx.enable(moderngl.DEPTH_TEST | moderngl.CULL_FACE | moderngl.BLEND)

        #gl2.glEnable(GL.GL_BLEND);
        #gl2.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        if not hasattr(self, 'programs'):
            #self.programs = [MainShader(), EdgeShader()]
            self.programs = [MainShader(), WireShader()]

        self.shaderkeys = {}
        for i, program in enumerate(self.programs):
            program.build(self.ctx).init()
            key = getattr(self.wnd.keys, f'F{i + 1}')
            self.shaderkeys[key] = self.toggleshader(program)

        self.materials = LazyMaterials(self.ctx, texture_directory)

        self.camera = Camera(self.wnd.keys, **kws)

        self.wnd.mouse_drag_event_func = self.camera.mouse_drag_event
        self.wnd.mouse_scroll_event_func = self.camera.mouse_scroll_event
        self.wnd.mouse_press_event_func = self.camera.mouse_press_event
        self.wnd.key_event_func = self.key_event

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
            elif key in self.shaderkeys:
                self.shaderkeys[key]()
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
        self.ctx.clear(*self.background)
        self.ctx.enable(moderngl.DEPTH_TEST)
        for program in self.programs:
            program.update(self.camera)
        for target in self.targets:
            target()

    def ___axes(self):
        # TODO: this defies the pattern of default shaders...
        # TODO: draw on top and with colors...
        axes = np.array([[[0, 0, 0], [1, 0, 0], [0, 0, 0]],
                         [[0, 0, 0], [0, 1, 0], [0, 0, 0]],
                         [[0, 0, 0], [0, 0, 1], [0, 0, 0]]])
        signature, program = self.programs[0]
        vbo = self.ctx.buffer(axes.astype('f4').tobytes())
        vao = self.ctx.simple_vertex_array(program, vbo, 'in_position')
        return (moderngl.LINES, None, vao)

    def draw_material(self, tf, material, T2F_N3F_V3F=None, wires=None):
        """Recursively aggregate data from a scene graph for a material"""
        T2F_N3F_V3F = [] if T2F_N3F_V3F is None else T2F_N3F_V3F
        wires = [] if wires is None else wires
        if hasattr(tf, 'models'):
            for model in tf.models:
                for mesh in model.meshes.get(material, ()):
                    try:
                        tridata = mesh.T2F_N3F_V3F(tf)
                        wiredata = mesh.wires(tf)
                        T2F_N3F_V3F.extend(tridata)
                        wires.extend(wiredata)
                    except:
                        print(f'failed to generate mesh {mesh}')
                        raise
        for child in tf.children:
            self.draw_material(child, material, T2F_N3F_V3F, wires)
        return T2F_N3F_V3F, wires

    def update(self):
        """Recompute the scenegraph and vao data (i.e. self.targets)"""
        if hasattr(self.source, 'scene'):
            self.scene = self.source.scene()
        else:
            self.scene = TForm()
            print(f'Source has no "scene" method: {self.source}')

        #self.targets = [self.axes()]

        self.targets = []
        for name in self.materials:
            vertexdata, wiredata = self.draw_material(self.scene, name)
            if vertexdata or wiredata:
                material = self.materials[name]
                for program in self.programs:
                    if program.active:
                        render = program.render(vertexdata, wiredata, material)
                        self.targets.append(render)

    @classmethod
    def test(cls, texture_directory, **kws):
        cube = Model.cube_model('generic_0', 1)
        source = Base(scene=(lambda : TForm(models=[cube])))
        inst = cls(source, texture_directory, **kws)
        inst.run()


def show(instance, texture_directory='../resources/textures', **kws):
    """Possible valid inputs:

        - Any object with a `scene` method
        - Tform instance
        - Mesh instance
        - Sequence of TForm instances
        - Sequence of Mesh instances

    """
    if not hasattr(instance, 'scene'):
        if isinstance(instance, TForm):
            tf = instance
        elif isinstance(instance, Mesh):
            model = Model(meshes={'generic_0': [instance]})
            tf = TForm(models=[model])
        elif all(isinstance(i, TForm) for i in instance):
            tf = TForm(children=instance)
        elif all(isinstance(i, Mesh) for i in instance):
            model = Model(meshes={'generic_0': instance})
            tf = TForm(models=[model])
        else:
            raise
        instance = Base(scene=lambda : tf)
    Window(instance, texture_directory, **kws).run()
