#pragma once
#include "basicparticle.h"

#include <vector>
#include <iostream>

class ParticleArchive {

public:
	ParticleArchive();
	ParticleArchive(const ParticleArchive& other);

	void addParticle(BasicParticle);
	// void addParticles(std::vector<MOOParticle>);
	void evaluateParticles();
	BasicParticle getParticle();
	std::vector<BasicParticle> getParticles();
	int size();

private:
	std::vector<BasicParticle> particles;

	void assignRanks();
	int dominates(BasicParticle& p1, BasicParticle& p2);
};