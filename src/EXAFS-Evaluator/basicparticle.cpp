#include "basicparticle.h"

BasicParticle::BasicParticle() {
	
}

BasicParticle::BasicParticle(std::vector< PDBAtom > atoms, double exafs_score, double potential_energy) {

	this->position = atoms;
	this->exafs_score = exafs_score;
	this->potential_energy = potential_energy;

	this->rank = -1;
}

BasicParticle::BasicParticle( const BasicParticle& other ) : position( other.position ),
															 exafs_score( other.exafs_score ),
															 potential_energy( other.potential_energy ),
															 rank (other.rank) { }