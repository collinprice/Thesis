#include "mooparticle.h"
#include "particlearchive.h"

#include <iostream>

MOOParticle::MOOParticle() : Chromosome() {

	// this->archive = new ParticleArchive();
}

MOOParticle::MOOParticle(double range, std::vector<PDBAtom> atoms) : Chromosome(atoms) {

	for (int i = 0; i < (int)atoms.size(); ++i) {
		velocity.push_back(PDBAtom(this->randomRange(range),this->randomRange(range),this->randomRange(range)));
	}

	// this->archive = new ParticleArchive();
}

// MOOParticle::~MOOParticle() {

// 	if (this->archive != NULL) {
// 		delete this->archive;
// 		this->archive = NULL;
// 	}
// }

MOOParticle::MOOParticle(const MOOParticle& other) 
: Chromosome(other),
  archive( other.archive ),
  velocity( other.velocity ) {

  	// ParticleArchive tempArchive = *other.archive;
  	// this->archive = new ParticleArchive(tempArchive);
}

void MOOParticle::updateBest() {
	// std::cout << "update pbest" << std::endl;
	this->archive.addParticle(this->getBasic());
	// std::cout << "eval pbest" << std::endl;
	this->archive.evaluateParticles();
}

void MOOParticle::updatePosition() {

	for (int i = 0; i < (int)this->atoms.size(); ++i) {
		this->atoms[i] = this->atoms[i] + this->velocity[i];
	}
}

void MOOParticle::updateVelocity(BasicParticle global_best, double inertia, double social, double cognitive) {
	
	std::vector<PDBAtom> personal_best = this->archive.getParticle().position;
	for (int i = 0; i < (int)this->atoms.size(); ++i) {

		this->velocity[i] = this->velocity[i]*inertia 
		+ (personal_best[i]-this->atoms[i])*social*this->unifRand()
		+ (global_best.position[i]-this->atoms[i])*cognitive*this->unifRand();
	}

	this->is_evaluated = false;
}

BasicParticle MOOParticle::getBasic() {
	return BasicParticle(this->atoms, this->exafs_score, this->potential_energy);
}

double MOOParticle::unifRand() {
	return rand() / double(RAND_MAX);
}

double MOOParticle::randomRange(double max) {

	double temp = max*this->unifRand();
	return this->unifRand() > 0.5 ? temp : -temp;
}