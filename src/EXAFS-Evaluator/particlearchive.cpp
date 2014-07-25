#include "particlearchive.h"
// #include "mooparticle.h"

ParticleArchive::ParticleArchive() {

}

ParticleArchive::ParticleArchive(const ParticleArchive& other) 
: particles( other.particles ) {}

void ParticleArchive::addParticle(BasicParticle particle) {
	// std::cout << "before addParticle" << std::endl;
	this->particles.push_back(particle);
	// std::cout << "after addParticle" << std::endl;
}

// void ParticleArchive::addParticles(std::vector<BasicParticle> particles) {

// 	for (std::vector<BasicParticle>::iterator iter = particles.begin(); iter != particles.end(); ++iter) {
// 		this->particles.push_back(*iter);
// 	}
// }

void ParticleArchive::evaluateParticles() {

	// std::cout << "Archive size before = " << this->particles.size() << std::endl;

	this->assignRanks();

	for (std::vector<BasicParticle>::iterator iter = this->particles.begin(); iter != this->particles.end(); /* skip ++iter */) {
		
		BasicParticle particle = *iter;
		// std::cout << "clean rank: " << particle.rank << std::endl;
		if (particle.rank != 1) {
			iter = this->particles.erase(iter);
		} else {
			++iter;
		}
	}

	// std::cout << "Archive size after = " << this->particles.size() << std::endl;
}

BasicParticle ParticleArchive::getParticle() {
	// std::cout << "Archive IN getParticle " << this->particles.size() << std::endl;	
	// Random index
	int index = (rand() / double(RAND_MAX)) * this->particles.size();
	// std::cout << "Archive getParticle at " << index << std::endl;	
	return this->particles.at(index);
}

std::vector<BasicParticle> ParticleArchive::getParticles() {
	return this->particles;
}

int ParticleArchive::size() {
	return this->particles.size();
}

void ParticleArchive::assignRanks() {

	// std::cout << "Assigning Ranks" << std::endl;

	// Catch if population size == 1
	if (this->particles.size() == 1) {
		this->particles[0].rank = 1;
		return;
	}

	int currentRank = 1;

	std::vector<BasicParticle> current_population = this->particles;
	std::vector<BasicParticle> new_population;

	// Reset all ranks.
	for (int i = 0; i < (int)current_population.size(); ++i) {
		current_population[i].rank = -1;
		// std::cout << current_population[i].exafs_score << " : " << current_population[i].potential_energy << std::endl;
	}
	
	while((int)current_population.size() > 0) {
		// std::cout << "SLOOP" << std::endl;
		int currentRankCount = 0;
		for (int i = 0; i < (int)current_population.size(); ++i) {
			
			bool should_keep = true;
			for (int j = 0; j < (int)current_population.size(); ++j) {
				
				if (i == j) continue;

				int result = this->dominates(current_population[i],current_population[j]);

				if (result == -1) {
					if (this->dominates(current_population[j],current_population[i]) == 1) {
						should_keep = false;
					}
				}
			}

			if (should_keep) {
				current_population[i].rank = currentRank;
				++currentRankCount;
			}
		}

		++currentRank;

		std::vector<BasicParticle> unranked_population;

		for (int i = 0; i < (int)current_population.size(); ++i) {
			// std::cout << "Rank: " << current_population[i].rank << std::endl;
			if (current_population[i].rank != -1) {
				new_population.push_back(current_population[i]);
			} else {
				unranked_population.push_back(current_population[i]);
			}
		}

		current_population = unranked_population;
		// std::cout << "ELOOP" << std::endl;
	}
	
	this->particles = new_population;
}

/*
	Checks is BasicParticle p1 dominate BasicParticle p2

	Returns:
		-1: no
		0: equal
		1: yes
*/
int ParticleArchive::dominates(BasicParticle& p1, BasicParticle& p2) {

	bool at_least_one_better = false;

	if (p1.exafs_score <= p2.exafs_score) {
		if (p1.exafs_score < p2.exafs_score) {
			at_least_one_better = true;
		}
	} else {
		return -1;
	}

	if (p1.potential_energy <= p2.potential_energy) {
		if (p1.potential_energy < p2.potential_energy) {
			at_least_one_better = true;
		}
	} else {
		return -1;
	}

	return at_least_one_better ? 1 : 0;
}