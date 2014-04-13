#ifndef HOUDINIFILEDUMPHELPER_H
#define HOUDINIFILEDUMPHELPER_H

#include <string>
#include <vector>
#include <ostream>

class HoudiniFileDumpParticles
{
public:
	class ParticlesDataProvider;

public:
	HoudiniFileDumpParticles();
	HoudiniFileDumpParticles(ParticlesDataProvider* dataProvider);
	~HoudiniFileDumpParticles();

	void setDataProvider(ParticlesDataProvider* dataProvider);
	void dump(std::ostream& stream);

private:
	ParticlesDataProvider*	_dataProvider;
};

class HoudiniFileDumpParticles::ParticlesDataProvider
{
public:

	ParticlesDataProvider() {}
	virtual ~ParticlesDataProvider() {}

	virtual int getNbPoints() = 0;
	virtual int getNbAttributes() = 0;
	virtual void getAttributesInfo(std::vector<std::string>& names, std::vector<int>& nbValues) = 0;

	virtual void getPtPosition(int ptID, float& posX, float& posY, float& posZ, float& posW) = 0;
	virtual void getPtAttributes(int ptID, const std::string& attribDelim, const std::string& valueDelim, std::ostream& out) = 0;
};

//------------------------------------------------------------------------------
// ConcreteDataProvider
//------------------------------------------------------------------------------
class houdini_Particle
{
public:
    float px, py, pz;
    float vx, vy, vz;
    float colorR, colorG, colorB;
    float mass;
};

//------------------------------------------------------------------------------
// ConcreteDataProvider
//------------------------------------------------------------------------------
class ConcreteDataProvider : public HoudiniFileDumpParticles::ParticlesDataProvider
    {
    private:
        houdini_Particle* _houdini_Particles;
        int	_nbhoudini_Particles;

    public:
        //-------------------------------------------------------------
        ConcreteDataProvider(houdini_Particle* houdini_Particles, int nbhoudini_Particles)
        : _houdini_Particles(houdini_Particles),
          _nbhoudini_Particles(nbhoudini_Particles)
        {}

        //-------------------------------------------------------------
    virtual ~ConcreteDataProvider()
        {}

        //-------------------------------------------------------------
        virtual int getNbPoints() { return _nbhoudini_Particles; }
        virtual int getNbAttributes() { return 3; }

        //-------------------------------------------------------------
        virtual void getAttributesInfo(std::vector<std::string>& names, std::vector<int>& nbValues)
        {
            // Spécifier le nom des attributs (tel qu'il apparaîtra dans houdini) et le nombre de valeur pour chaque attribut.
            // Un simple float nécessitera une seule valeur, alors que vecteur (par exemple la vitesse) nécessitera 3 valeurs.
            // Note: La position n'est pas considérée comme un attribut!
            names.push_back("v");
            nbValues.push_back(3);	// v = (vx, vy, vz)

            names.push_back("color");
            nbValues.push_back(3);	// color = (colorR, colorG, colorB)

            names.push_back("mass");
            nbValues.push_back(1); // mass = (mass)
        }

        //-------------------------------------------------------------
        virtual void getPtPosition(int ptID, float& posX, float& posY, float& posZ, float& posW)
        {
            posX = _houdini_Particles[ptID].px;
            posY = _houdini_Particles[ptID].py;
            posZ = _houdini_Particles[ptID].pz;
            posW = 0; //Not set for now
        }

        //-------------------------------------------------------------
        virtual void getPtAttributes(int ptID, const std::string& attribDelim, const std::string& valueDelim, std::ostream& out)
        {
            // Dans cette fonction, la valeur des attributs est outputé. Ceci permet de ne pas être dépendant des types utilisés pour
            // stocker les attributs.

            // attribDelim et valueDelim représentent respectivement le/les caractère(s) d'espacement entre deux attributs
            // et deux valeurs d'un même attribut.

            // Note: Les attributs doivent être spécifiés dans le même ordre que dans la fonction getAttributesInfo
            out << _houdini_Particles[ptID].vx << valueDelim << _houdini_Particles[ptID].vy << valueDelim << _houdini_Particles[ptID].vz << attribDelim;
            out << _houdini_Particles[ptID].colorR << valueDelim << _houdini_Particles[ptID].colorG << valueDelim << _houdini_Particles[ptID].colorB << attribDelim;
            out << _houdini_Particles[ptID].mass; // Pas besoin de délimiteur car c'est le dernier attribut!
        }
    };

#endif	// HOUDINIFILEDUMPHELPER_H
