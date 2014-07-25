#include "moopso.h"

#include <iostream>
#include <math.h>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>

MOOPSO::MOOPSO(EXAFSEvaluator* exafs_evaluator, double inertia, double social, double cognitive, double velocity_range, int max_generations, std::string results_file) {

	this->exafs_evaluator = exafs_evaluator;
	this->velocity_range = velocity_range;
	this->inertia = inertia;
	this->social = social;
	this->cognitive = cognitive;
	this->max_generations = max_generations;
	this->results_file = results_file;

	// this->archive = new ParticleArchive();
}

MOOPSO::~MOOPSO() {

	delete this->exafs_evaluator;
	this->exafs_evaluator = NULL;

	// delete this->archive;
	// this->archive = NULL;
}

bool chromosome_sort_moopso(MOOParticle const & a, MOOParticle const & b) {
	return a.rank < b.rank;
}

void MOOPSO::begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations) {

	for (int i = 0; i < (int)initial_populations.size(); ++i) {
		
		std::stringstream ss;
		ss << (i+1);
		this->stats_folder = "run" + ss.str();
		mkdir(this->stats_folder.c_str(), 0755);

		this->initPopulation(initial_populations[i]);
		this->initStats();
		this->recordStats();

		// this->global_best_particle = this->archive->getParticle();
		this->updateGlobalBest();

		std::cout << "Begin Run " << (i+1) << std::endl;
		for (int i = 0; i < this->max_generations; ++i) {
			std::cout << "Generation: " << (i+1) << std::endl;

			this->updateVelocities();
			this->updatePositions();

			this->evaluatePopulation();
			this->recordStats();

			this->updateGlobalBest();

			if (this->convergence()) break;
		}

		this->finalStats();
	}
}

void MOOPSO::updateVelocities() {

	for (int i = 0; i < (int)this->population.size(); ++i) {
		this->population[i].updateVelocity(this->global_best_particle,this->inertia,this->social,this->cognitive);
	}
}

void MOOPSO::updatePositions() {

	for (int i = 0; i < (int)this->population.size(); ++i) {
		this->population[i].updatePosition();
	}
}

void MOOPSO::updateGlobalBest() {

	// std::cout << "updateGlobalBest" << std::endl;
	for (int i = 0; i < (int)this->population.size(); ++i) {

		// std::cout << "befor wtf" << std::endl;
		this->archive.addParticle(this->population[i].getBasic());
		// std::cout << "after wtf" << std::endl;

		this->population[i].updateBest();
		// std::cout << "after updateBest" << std::endl;
		// std::cout << "after updateBest" << std::endl;
		// std::cout << "after updateBest" << std::endl;
		// std::cout << "after updateBest" << std::endl;
		// std::cout << "after updateBest" << std::endl;
		// std::cout << "after updateBest" << std::endl;

		// this->archive->addParticle(this->population[i]);
		// std::cout << "after addParticle to global archive" << std::endl;
	}

	// std::cout << "before eval gbest" << std::endl;
	this->archive.evaluateParticles();
	// std::cout << "after eval gbest" << std::endl;
	this->global_best_particle = this->archive.getParticle();

	// std::cout << "Global Archive size = " << this->archive.size() << std::endl;
}

void MOOPSO::initPopulation(std::vector< std::vector<PDBAtom> > population) {

	// std::cout << "initPopulation" << std::endl;
	this->population.clear();
	for (std::vector< std::vector<PDBAtom> >::iterator i = population.begin(); i != population.end(); ++i) {
		
		MOOParticle child(this->velocity_range,*i);
		this->population.push_back(child);
	}

	// std::cout << "end initPopulation" << std::endl;

	this->evaluatePopulation();

	// std::cout << "Archive Test: Start" << std::endl;

	// this->archive.addParticle(this->population[0].getBasic());
	// this->archive.addParticle(this->population[1].getBasic());
	// this->archive.addParticle(this->population[2].getBasic());
	// this->archive.addParticle(this->population[3].getBasic());
	// this->archive.addParticle(this->population[4].getBasic());

	// this->archive.evaluateParticles();

	// std::cout << "Archive Test: Finish" << std::endl;

	this->updateGlobalBest();
}

void MOOPSO::evaluatePopulation() {
	for (int i = 0; i < (int)this->population.size(); ++i) {
		if (!this->population[i].is_evaluated) {
			this->evaluate(this->population[i]);
			this->population[i].is_evaluated = true;
		}
		std::cout << "\t Child: " << i << ", " << this->population[i].exafs_score << ", " << this->population[i].potential_energy << std::endl;
	}
}

void MOOPSO::evaluate( MOOParticle& child ) {

	this->exafs_evaluator->updateAtoms(child.atoms);
	
	child.exafs_score = this->exafs_evaluator->calculateRMSD();
	child.exafs_data = this->exafs_evaluator->getEXAFSData();

	child.potential_energy = this->exafs_evaluator->calculatePotentialEnergy();
}

// Particle MOOPSO::best_particle() {

// 	double best = std::numeric_limits<double>::max();
// 	Particle best_particle;
// 	for (std::vector<Particle>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
// 		if (child->exafs_score < best) {
// 			best = child->exafs_score;
// 			best_particle = *child;
// 		}
// 	}

// 	return best_particle;
// }

// void MOOPSO::saveBestParticle() {

// 	Particle best_individual = this->best_particle();
// 	this->best_individuals.push_back(best_individual);

// 	std::cout << "Best: " << best_individual.exafs_score << std::endl;
// 	std::cout << "Global Best: " << global_best_particle.exafs_score << std::endl;
// }

double MOOPSO::unifRand() {
	return rand() / double(RAND_MAX);
}

double MOOPSO::randInt(int max) {

	return this->unifRand() * max;
}

void MOOPSO::initStats() {

	this->output_stream.open(( this->stats_folder + "/" + this->results_file).c_str());
	if (this->output_stream.is_open()) {
		std::cout << "Results file ready." << std::endl;
	} else {
		std::cout << "Could not create results file." << std::endl;
	}
}

void MOOPSO::recordStats() {

	std::vector<BasicParticle> best_particles = this->archive.getParticles();

	std::stringstream exafs_score_string;
	std::stringstream potential_energy_string;

	for (std::vector<BasicParticle>::iterator iter = best_particles.begin(); iter != best_particles.end(); ++iter) {
		
		exafs_score_string << iter->exafs_score << ",";
		potential_energy_string << iter->potential_energy << ",";
	}

	this->output_stream << exafs_score_string.str() << std::endl;
	this->output_stream << potential_energy_string.str() << std::endl;
	// double average_exafs_score = 0;
	// double best_score = std::numeric_limits<double>::max();
	// Particle best_chromosome;

	// for (std::vector<Particle>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
	// 	average_exafs_score += child->exafs_score;
	// 	if (child->exafs_score < best_score) {
	// 		best_score = child->exafs_score;
	// 		best_chromosome = *child;
	// 	}
	// }

	// this->output_stream << best_score << "," << (average_exafs_score/(int)this->population.size()) << "," << best_chromosome.potential_energy << std::endl;
}

void MOOPSO::finalStats() {

	std::cout << "Final Stats." << std::endl;

	this->output_stream.close();

	// TODO: Bring back later!!!
	// std::stringstream output_chromosome_filename;
	// int counter = 0;
	// std::vector<MOOParticle> best_particles = this->archive->getParticles();
	// for (std::vector<MOOParticle>::iterator child = best_particles.begin(); child != best_particles.end(); ++child) {
		
	// 	output_chromosome_filename << "/best-chromosome-" << (counter++) << ".pdb";
	// 	std::cout << "Writing out: " << this->stats_folder << output_chromosome_filename.str() << std::endl;
		
	// 	this->exafs_evaluator->updateAtoms(child->atoms);
	// 	this->exafs_evaluator->writePDB(this->stats_folder + output_chromosome_filename.str());
	// 	output_chromosome_filename.str("");
	// }

	// this->output_stream.close();

	// std::vector< std::pair<double, double> > target_exafs = this->exafs_evaluator->getTargetEXAFS();
	// std::ofstream output(this->stats_folder + "/generation_data.csv");
	// for (int i = 0; i < (int)this->best_individuals[0].exafs_data.size(); ++i) {
		
	// 	output << this->best_individuals[0].exafs_data[i].first;

	// 	if (this->best_individuals[0].exafs_data[i].first != 0) {
	// 		output << "," << target_exafs[i-1].second;
	// 	} else {
	// 		output << ",0";
	// 	}
		
	// 	for (int j = 0; j < (int)this->best_individuals.size(); ++j) {
			
	// 		if ((int)this->best_individuals[j].exafs_data.size() == j) continue; // guard for no exafs data. occurs when no individuals are valid.

	// 		output << "," << this->best_individuals[j].exafs_data[i].second;
	// 	}
	// 	output << std::endl;
	// }
	// output.close();

	// this->exafs_evaluator->updateAtoms(this->best_individuals[this->best_individuals.size()-1].atoms);
	// this->exafs_evaluator->writePDB(this->stats_folder + "/best_chromosome.pdb");
	
}

bool MOOPSO::convergence() {

	std::vector<MOOParticle> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort_moopso);
	std::vector<MOOParticle>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) == 1;
}

bool MOOPSO::convergence(double rate) {

	std::vector<MOOParticle> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort_moopso);
	std::vector<MOOParticle>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) <= (rate * this->population.size());
}
