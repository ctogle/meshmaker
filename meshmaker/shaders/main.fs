#version 330

struct Material {
    sampler2D diffuse;
    vec3      specular;
    float     shininess;
};


struct DirLight {
    vec3 direction;
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};


uniform Material material;
uniform DirLight dirLight1;
uniform DirLight dirLight2;
uniform vec3 viewDir;


in vec3 v_vert;
in vec3 v_norm;
in vec2 texcoord_0;


out vec4 f_color;


vec3 CalcDirLight(DirLight light, vec3 normal, vec3 viewDir)
{
    vec3 direction = normalize(-light.direction);
    float lum = clamp(dot(direction, normal), 0.1, 1.0);
    float spec = pow(max(dot(viewDir, reflect(direction, normal)), 0.02), material.shininess);

    vec3 difftex = vec3(texture(material.diffuse, texcoord_0));
    //vec3 spectex = vec3(texture(material.specular, texcoord_0));
    vec3 spectex = material.specular;
    
    vec3 ambient = light.ambient * difftex;
    vec3 diffuse = light.diffuse * lum * difftex;
    vec3 specular = light.specular * spec * spectex;    
    return (ambient + diffuse + specular - specular);
}


void main() {
    f_color = vec4(0.0);
    f_color += vec4(CalcDirLight(dirLight1, normalize(v_norm), normalize(viewDir)), 1.0);
    f_color += vec4(CalcDirLight(dirLight2, normalize(v_norm), normalize(viewDir)), 1.0);
    f_color.a = 1.0;
}
