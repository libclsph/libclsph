#include "Partio.h"

//Code is inspired by the examples provided with partio
//We try to not change the style of the code to not break anything
inline Partio::ParticlesDataMutable* makeData(particle* buffer,simulation_parameters parameters) {
 
    Partio::ParticlesDataMutable& data=*Partio::create();
    Partio::ParticleAttribute positionAttr=data.addAttribute("position",Partio::VECTOR,3);
    Partio::ParticleAttribute velocityAttr=data.addAttribute("velocity",Partio::VECTOR,3);
    Partio::ParticleAttribute colorAttr=data.addAttribute("color",Partio::VECTOR,3);
    Partio::ParticleAttribute idAttr=data.addAttribute("id",Partio::INT,1);
    Partio::ParticleAttribute massAttr = data.addAttribute("mass", Partio::FLOAT, 1);
 
    for(int i=0; i<parameters.particles_count; ++i) {
 
        particle temp_particle = buffer[i];
 
        Partio::ParticleIndex index=data.addParticle();
 
        float* pos=data.dataWrite<float>(positionAttr,index);
        float* vel=data.dataWrite<float>(velocityAttr,index);
        float* color=data.dataWrite<float>(colorAttr,index);
 
        float* mass=data.dataWrite<float>(massAttr,index);
 
        int* id=data.dataWrite<int>(idAttr,index);
 
        pos[0]= temp_particle.position.s[0];
        pos[1]= temp_particle.position.s[1];
        pos[2]= temp_particle.position.s[2];
 
        vel[0]= temp_particle.velocity.s[0];
        vel[0]= temp_particle.velocity.s[0];
        vel[0]= temp_particle.velocity.s[0];
 
        //TODO Use actual colors---------------------------------------------------------
        color[0] =
            temp_particle.density > 1000.f && temp_particle.density <= 2000.f ?
            (temp_particle.density - 1000.f) / 1000.f :
            0.f;
 
        color[1] =
            temp_particle.density >= 0.f && temp_particle.density < 1000.f ?
            1.f - (temp_particle.density) / 1000.f :
            0.f;
 
        color[2] =
            temp_particle.density >= 500.f && temp_particle.density <= 1000.f ?
            (temp_particle.density - 500.f) / 500.f :
            temp_particle.density >= 1000.f && temp_particle.density <= 1500.f ?
            1.f - (temp_particle.density - 1000.f) / 500.f:
            0.f;
        //-------------------------------------------------------------------------------
 
        mass[0]=parameters.particle_mass;
 
        id[0]=index;
 
    }
    return &data;
}
