#include "HoudiniFileDumpHelper.h"

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// HoudiniFileDumpParticles ++++++++++++++++++++++++++++++++++++++++++++++++++++

// -----------------------------------------------------------------------------
// Constructors / Destructor
HoudiniFileDumpParticles::HoudiniFileDumpParticles()
: _dataProvider(0x0)
{}

HoudiniFileDumpParticles::HoudiniFileDumpParticles(ParticlesDataProvider* dataProvider)
: _dataProvider(dataProvider)
{}

HoudiniFileDumpParticles::~HoudiniFileDumpParticles()
{}

void HoudiniFileDumpParticles::dump(std::ostream& stream)
{
	if (!_dataProvider) { return; }

	// Write header
	int nbPoints = _dataProvider->getNbPoints();
	int nbAttributes = _dataProvider->getNbAttributes();
	stream << "PGEOMETRY V5" << std::endl;
	stream << "NPoints " << _dataProvider->getNbPoints() << " NPrims 1" << std::endl;
	stream << "NPointGroups 0 NPrimGroups 1" << std::endl;
	stream << "NPointAttrib " << nbAttributes << " NVertexAttrib 0 NPrimAttrib 2 NAttrib 0" << std::endl;
	
	// Attribs
	std::vector<std::string>	attribNames;
	std::vector<int>			attribNbValues;
	_dataProvider->getAttributesInfo(attribNames, attribNbValues);
	stream << "PointAttrib" << std::endl;
	for (int i=0; i<nbAttributes; ++i)
	{
		stream << attribNames[i] << " " << attribNbValues[i] << " float";	// TODO: Supporter d'autres types que les floats!!!'
		for (int v=0; v<attribNbValues[i]; ++v)
		{
			stream << " 1";
		}
		stream << std::endl;
	}

	// Dump points
	float px, py, pz, pw;
	std::vector<float> values;
	for (int n=0; n<nbPoints; ++n)
	{
		// Get Pt Datas
		values.clear();
		_dataProvider->getPtPosition(n, px, py, pz, pw);

		// Output point position
		stream << px << " " << py << " " << pz << " " << pw << " ";

		// Output point attributes
		if (nbAttributes > 0)
		{
			stream << "(";
			_dataProvider->getPtAttributes(n, "\t", " ", stream);
			stream << ")" << std::endl;
		}		
	}

	// Primitive Attribs
	stream << "PrimitiveAttrib" << std::endl;
	stream << "generator 1 index 1 location1" << std::endl;
	stream << "dopobject 1 index 1 /obj/AutoDopNetwork:1" << std::endl;

	// Particles
	stream << "Part " << nbPoints;
	for (int n=0; n<nbPoints; ++n)
	{
		stream << " " << n;
	}
	stream << " [0\t0]" << std::endl;

	// Detail attributes
	// Nada!

	// ???
	stream << "box_object1 unordered" << std::endl;
	stream << "1 1" << std::endl;

	// Extra
	stream << "beginExtra" << std::endl;
	stream << "endExtra" << std::endl;
}

// -----------------------------------------------------------------------------
// Public functions
void HoudiniFileDumpParticles::setDataProvider(ParticlesDataProvider* dataProvider)
{
	_dataProvider = dataProvider;
}