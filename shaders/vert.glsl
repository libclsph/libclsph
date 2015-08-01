#version 330 core

out vec3 color;
in vec3 position;
in float density;

uniform float base_density;
uniform mat4 model_view_projection;

void main() {
  float minimum = base_density * 0.5f;
  float maximum = base_density * 1.5f;

  color = vec3(1.f - ((density - base_density) / base_density),
               ((density - base_density) / base_density), 0.f);

   gl_Position = model_view_projection * vec4(position, 1.0);
}
