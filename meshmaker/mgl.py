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

    def __init__(self,
                 vs='../meshmaker/shaders/main.vs', gs=None,
                 fs='../meshmaker/shaders/main.fs',
                 signature=('in_texcoord_0', 'in_normal', 'in_position')):
        self.vs = vs
        self.gs = gs
        self.fs = fs
        self.signature = signature

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
        self._program["dirLight.direction"].value = (-2.0,-1.0,-3.0)
        self._program["dirLight.diffuse"].value   = ( 1.0, 1.0, 1.0)
        self._program["dirLight.ambient"].value   = ( 0.2, 0.2, 0.2)
        self._program["dirLight.specular"].value  = ( 0.8, 0.8, 0.8)
        self._program["material.specular"].value  = ( 0.9, 0.9, 0.9)
        self._program["material.shininess"].value = 100

    def update(self, camera):
        """Set GLSL uniform variables per frame for camera"""
        self._program['projection'].write(camera.projection.astype('f4').tobytes())
        self._program['view'].write(camera.view.astype('f4').tobytes())
        self._program['viewDir'].value = tuple(camera.viewdir)

    def render(self, vertexdata, material):
        """Compute function which renders vertexdata using shaders"""
        vertexdata = np.array(vertexdata)
        vbo = self._ctx.buffer(vertexdata.astype('f4').tobytes())
        vao = self._ctx.simple_vertex_array(self._program, vbo, *self.signature)

        def f():
            if material:
                material.use()
            vao.render(moderngl.TRIANGLES)

        return f


class Window(Base):

    defaultshaders = {
        'signature': ('in_texcoord_0', 'in_normal', 'in_position'),
        'vertex_shader': '''
            #version 330

            uniform mat4 projection;
            uniform mat4 view;

            in vec3 in_position;
            in vec3 in_normal;
            in vec2 in_texcoord_0;

            out vec3 v_vert;
            out vec3 v_norm;
            out vec2 v_text;

            void main() {
                gl_Position = (projection * view) * vec4(in_position, 1.0);
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
        super().__init__(**kws)
        self.source = source
        self.targets = []

        settings.WINDOW['class'] = 'moderngl_window.context.pyglet.Window'
        self.wnd = moderngl_window.create_window_from_settings()

        self.ctx = self.wnd.ctx
        self.ctx.enable(moderngl.DEPTH_TEST | moderngl.CULL_FACE)

        #if not hasattr(self, 'shaders'):
        #    self.shaders = [self.defaultshaders.copy()]

        if not hasattr(self, 'programs'):
            vs = os.path.join(os.path.dirname(__file__), 'shaders', 'main.vs')
            fs = os.path.join(os.path.dirname(__file__), 'shaders', 'main.fs')
            self.programs = [ShaderProgram(vs=vs, fs=fs)]

        for program in self.programs:
            program.build(self.ctx).init()

        '''
        self.programs = []
        for spec in self.shaders:
            callback = spec.pop('callback', None)
            signature = spec.pop('signature',
                ('in_texcoord_0', 'in_normal', 'in_position'))
            program = self.ctx.program(**spec)
            if callback:
                callback(program)
            self.programs.append((signature, program))
        '''

        '''
        callback = self.shaders.pop('callback', None)
        signature = self.shaders.pop('signature',
            ('in_texcoord_0', 'in_normal', 'in_position'))
        self.programs = [(signature, self.ctx.program(**self.shaders))]
        if callback:
            for signature, program in self.programs:
                callback(program)
        '''


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

        for program in self.programs:
            program.update(self.camera)

        #for signature, program in self.programs:
        #    program['projection'].write(self.camera.projection.astype('f4').tobytes())
        #    program['view'].write(self.camera.view.astype('f4').tobytes())

        for target in self.targets:
            target()

        #for mode, texture, vao in self.targets:
        #    if texture is not None:
        #        texture.use()
        #    vao.render(mode)

    def axes(self):
        # TODO: this defies the pattern of default shaders...
        # TODO: draw on top and with colors...
        axes = np.array([[[0, 0, 0], [1, 0, 0], [0, 0, 0]],
                         [[0, 0, 0], [0, 1, 0], [0, 0, 0]],
                         [[0, 0, 0], [0, 0, 1], [0, 0, 0]]])
        signature, program = self.programs[0]
        vbo = self.ctx.buffer(axes.astype('f4').tobytes())
        vao = self.ctx.simple_vertex_array(program, vbo, 'in_position')
        return (moderngl.LINES, None, vao)

    def draw_material(self, tf, material, T2F_N3F_V3F=None, wireframe=None):
        """Recursively aggregate data from a scene graph for a material"""
        T2F_N3F_V3F = [] if T2F_N3F_V3F is None else T2F_N3F_V3F
        wireframe = [] if wireframe is None else wireframe
        if hasattr(tf, 'models'):
            for model in tf.models:
                for mesh in model.meshes.get(material, ()):
                    try:
                        tridata, wiredata = mesh.T2F_N3F_V3F(tf)
                        T2F_N3F_V3F.extend(tridata)
                        wireframe.extend(wiredata)
                    except:
                        print(f'failed to generate mesh {mesh}')
                        raise
        for child in tf.children:
            self.draw_material(child, material, T2F_N3F_V3F, wireframe)
        return T2F_N3F_V3F, wireframe

    def update(self):
        """Recompute the scenegraph and vao data (i.e. self.targets)"""
        if hasattr(self.source, 'scene'):
            self.scene = self.source.scene()
        else:
            self.scene = TForm()
            print(f'Source has no "scene" method: {self.source}')

        #self.targets = [self.axes()]
        self.targets = []

        # vertex_data -> vbo
        # program + vbo -> vao
        # mode + material + vao -> material.use;vao.render(mode)

        # mesh could generate vao,mode,material from program
        # e.g. T2F_N3F_V3F specific shader program

        # must compute vertex data on a per material basis
        # can reuse vertex data between shader programs via adaptor

        for name in self.materials:
            vertexdata, _ = self.draw_material(self.scene, name)
            if vertexdata:
                material = self.materials[name]
                for program in self.programs:
                    render = program.render(vertexdata, material)
                    self.targets.append(render)

        '''
        for signature, program in self.programs:
            for name in self.materials:
                vertices, wireframe = self.draw_material(self.scene, name)
                if vertices:
                    if len(signature) == 2:
                        mat = None
                        notuv = lambda i: not ((i % 8 == 0) or (i % 8 == 1))
                        vertices = [v for i, v in enumerate(vertices) if notuv(i)]
                    elif len(signature) == 3:
                        mat = self.materials[name]
                    vertices = np.array(vertices)
                    vbo = self.ctx.buffer(vertices.astype('f4').tobytes())
                    vao = self.ctx.simple_vertex_array(program, vbo, *signature)
                    self.targets.append((moderngl.TRIANGLES, mat, vao))
                    #'in_texcoord_0', 'in_normal', 'in_position')
                    #'in_position', 'in_normal')
                    #'in_position', 'in_normal', 'in_texcoord_0')
        '''

    @classmethod
    def test(cls, texture_directory, **kws):
        cube = Model.cube_model('generic_0', 1)
        source = Base(scene=(lambda : TForm(models=[cube])))
        inst = cls(source, texture_directory, **kws)
        inst.run()


def show(instance, texture_directory='../resources/textures', **kws):
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
