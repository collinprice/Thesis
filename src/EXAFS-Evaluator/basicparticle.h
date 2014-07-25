#pragma once
#include "pdbatom.h"

#include <vector>

class BasicParticle {

public:
	BasicParticle();
	BasicParticle(std::vector< PDBAtom >, double, double);
	BasicParticle( const BasicParticle& other );

	std::vector< PDBAtom > position;

	double exafs_score;
	double potential_energy;
	int rank;
	
};