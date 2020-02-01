#version 330

uniform vec3 viewDir;

uniform sampler2D Texture;

struct Material {
    sampler2D diffuse;
    vec3      specular;
    float     shininess;
};
uniform Material material;

struct DirLight {
    vec3 direction;
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};
uniform DirLight dirLight;

in vec3 v_vert;
in vec3 v_norm;
in vec2 texcoord_0;

out vec4 f_color;


vec3 CalcDirLight(DirLight light, vec3 normal, vec3 viewDir)
{
    vec3 direction = normalize(-light.direction);
    float lum = clamp(dot(direction, normal), 0.0, 1.0);
    float spec = pow(max(dot(viewDir, reflect(direction, normal)), 0.0), material.shininess);

    vec3 difftex = vec3(texture(material.diffuse, texcoord_0));
    //vec3 spectex = vec3(texture(material.specular, texcoord_0));
    vec3 spectex = material.specular;
    
    vec3 ambient = light.ambient * difftex;
    vec3 diffuse = light.diffuse * lum * difftex;
    vec3 specular = light.specular * spec * spectex;    
    return (ambient + diffuse + specular);
}


void main() {
    f_color = vec4(CalcDirLight(dirLight, normalize(v_norm), normalize(viewDir)), 1.0);
}
