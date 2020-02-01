#version 330

uniform mat4 projection;
uniform mat4 view;

in vec3 in_position;
in vec3 in_normal;
in vec2 in_texcoord_0;

out vec3 v_vert;
out vec3 v_norm;
out vec2 texcoord_0;

void main() {
    gl_Position = (projection * view) * vec4(in_position, 1.0);
    v_vert = in_position;
    v_norm = in_normal;
    texcoord_0 = in_texcoord_0;
}
