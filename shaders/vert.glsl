#version 330 core

out vec3 color;
in vec3 position;
in float density;

uniform float base_density;

void main() {
    
   if(density > base_density){
    color=vec3(0.0,1.0,0.0);
   }
   else{
    color=vec3(1.0,0.0,0.0);
   }
   
   
   gl_Position = vec4(position, 1.0);
}
