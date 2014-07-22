#include "exafsmoode.h"

#include <iostream>
#include <math.h>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>

EXAFSMOODE::EXAFSMOODE(EXAFSEvaluator* exafs_evaluator, double f, double cr, int max_generations, std::string results_file) {

	this->exafs_evaluator = exafs_evaluator;
	this->f = f;
	this->cr = cr;
	this->max_generations = max_generations;
	this->results_file = results_file;
}

EXAFSMOODE::~EXAFSMOODE() {

	delete this->exafs_evaluator;
}

bool chromosome_sort_moode(Chromosome const & a, Chromosome const & b) {
	return a.rank < b.rank;
}

void EXAFSMOODE::begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations) {

	for (int i = 0; i < (int)initial_populations.size(); ++i) {
		
		std::stringstream ss;
		ss << (i+1);
		this->stats_folder = "run" + ss.str();
		mkdir(this->stats_folder.c_str(), 0755);

		this->initPopulation(initial_populations[i]);
		this->initStats();
		this->recordStats();

		std::cout << "Begin Run " << (i+1) << std::endl;
		for (int i = 0; i < this->max_generations; ++i) {

			std::cout << "Generation: " << (i+1) << std::endl;
			this->evolve();

			if (this->convergence()) break;
		}

		this->finalStats();
	}
}

void EXAFSMOODE::begin_recentering(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations, int population_size, double convergence_rate, int max_iterations) {

	for (int i = 0; i < (int)initial_populations.size(); ++i) {

		std::cout << "Begin Run " << (i+1) << std::endl;

		std::stringstream ss;
		ss << (i+1);
		this->stats_folder = "run" + ss.str();
		mkdir(this->stats_folder.c_str(), 0755);

		// Get subset of population.
		std::random_shuffle(initial_populations[i].begin(), initial_populations[i].end());
		std::vector< std::vector<PDBAtom> > subset_pop = std::vector< std::vector<PDBAtom> >(initial_populations[i].begin(), initial_populations[i].begin()+population_size);

		this->initPopulation(subset_pop);
		this->initStats();
		this->recordStats();

		for (int j = 0; j < max_iterations; ++j) {
			
			std::cout << "Iteration: " << (j+1) << std::endl;
			for (int epoch = 0; epoch < this->max_generations; ++epoch) {
				
				std::cout << "Generation: " << (epoch+1) << std::endl;
				this->evolve();

				if (j != (max_iterations-1) && this->convergence(convergence_rate)) {

					std::cout << "Converged" << std::endl;

					std::vector<Chromosome> population_copy = this->population;
					std::sort(population_copy.begin(), population_copy.end(), chromosome_sort_moode);
					std::vector<Chromosome>::iterator it = std::unique(population_copy.begin(), population_copy.end());
					population_copy.resize( std::distance(population_copy.begin(), it) );

					// Shuffle library again and get subset minus best individuals
					std::random_shuffle(initial_populations[i].begin(), initial_populations[i].end());
					subset_pop = std::vector< std::vector<PDBAtom> >(initial_populations[i].begin(), initial_populations[i].begin()+(population_size-(int)population_copy.size()));

					for (std::vector<Chromosome>::iterator iter = population_copy.begin(); iter != population_copy.end(); ++iter) {
						subset_pop.push_back(iter->atoms);
					}

					this->initPopulation(subset_pop);
					this->recordStats();
					break;
				}
			}
		}

		this->finalStats();
	}
}

void EXAFSMOODE::initPopulation(std::vector< std::vector<PDBAtom> > population) {

	this->best_individuals.clear();
	this->population.clear();
	for (std::vector< std::vector<PDBAtom> >::iterator i = population.begin(); i != population.end(); ++i) {
		
		Chromosome child(*i);
		this->population.push_back(child);
	}

	this->evaluatePopulation();
}

void EXAFSMOODE::evaluatePopulation() {
	for (int i = 0; i < (int)this->population.size(); ++i) {
		if (!this->population[i].is_evaluated) {
			this->evaluate(this->population[i]);
			this->population[i].is_evaluated = true;
		}
		std::cout << "\t Child: " << i << ", " << this->population[i].exafs_score << ", " << this->population[i].potential_energy << std::endl;
	}

	this->rankPopulation();
}

void EXAFSMOODE::evaluate( Chromosome& child ) {

	this->exafs_evaluator->updateAtoms(child.atoms);
	
	child.exafs_score = this->exafs_evaluator->calculateRMSD();
	child.exafs_data = this->exafs_evaluator->getEXAFSData();

	child.potential_energy = this->exafs_evaluator->calculatePotentialEnergy();
}

void EXAFSMOODE::evolve() {
	
	std::vector<Chromosome> new_population;

	for (int i = 0; i < (int)this->population.size(); ++i) {
		
		int r1;
		int r2;
		int r3;

		do {
			r1 = this->randInt((int)this->population.size());
		} while (i == r1);

		do {
			r2 = this->randInt((int)this->population.size());
		} while (i == r2 || r1 == r2);

		do {
			r3 = this->randInt((int)this->population.size());
		} while (i == r3 || r1 == r3 || r2 == r3);

		new_population.push_back(this->mutate(i,r1,r2,r3));
	}

	this->population = new_population;
	this->rankPopulation();
	this->recordStats();
}

Chromosome EXAFSMOODE::mutate(int i, int r1, int r2, int r3) {

	std::vector<PDBAtom> modified_atoms;
	Chromosome selectedChromosome = this->population[i];

	int chromosome_length = selectedChromosome.atoms.size();

	// Modify each atom.
	for (int j = 0; j < chromosome_length; ++j) {
		
		PDBAtom temp = this->population[r2].atoms[j] - this->population[r3].atoms[j];
		temp = temp * this->f;
		modified_atoms.push_back(this->population[r1].atoms[j] + temp);
	}

	int Irand = this->randInt(chromosome_length);
	std::vector<PDBAtom> final_atoms;

	// Decide which atoms to replace.
	for (int j = 0; j < chromosome_length; ++j) {
		
		if (this->unifRand() <= this->cr || j == Irand) {
			final_atoms.push_back(modified_atoms[j]);
		} else {
			final_atoms.push_back(selectedChromosome.atoms[j]);
		}
	}

	// Test if fitness has improved.
	Chromosome originalChromosome = selectedChromosome;

	// Replace atoms and evaluate.
	selectedChromosome.atoms = final_atoms;
	this->evaluate(selectedChromosome);

	if (this->chromosomeDominates(selectedChromosome, originalChromosome) > -1) {
		std::cout << "\t Child: " << i << ", Modified: " << selectedChromosome.exafs_score << ", " << selectedChromosome.potential_energy << std::endl;
		return selectedChromosome;
	} else {
		std::cout << "\t Child: " << i << ", Original: " << originalChromosome.exafs_score << ", " << originalChromosome.potential_energy << std::endl;
		return originalChromosome;
	}
}

Chromosome EXAFSMOODE::best_chromosome() {

	double best = std::numeric_limits<double>::max();
	Chromosome best_chromosome;
	for (std::vector<Chromosome>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		if (child->exafs_score < best) {
			best = child->exafs_score;
			best_chromosome = *child;
		}
	}

	return best_chromosome;
}

void EXAFSMOODE::saveBestChromosome() {

	Chromosome best_individual = this->best_chromosome();
	this->best_individuals.push_back(best_individual);

	std::cout << "Best: " << best_individual.exafs_score << std::endl;
}

double EXAFSMOODE::unifRand() {
	return rand() / double(RAND_MAX);
}

double EXAFSMOODE::randInt(int max) {

	return this->unifRand() * max;
}

void EXAFSMOODE::initStats() {

	this->output_stream.open(( this->stats_folder + "/" + this->results_file).c_str());
	if (this->output_stream.is_open()) {
		std::cout << "Results file ready." << std::endl;
	} else {
		std::cout << "Could not create results file." << std::endl;
	}
}

void EXAFSMOODE::recordStats() {

	std::stringstream exafs_score_string;
	std::stringstream potential_energy_string;

	for (std::vector<Chromosome>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		if (child->rank == 1) {
			
			exafs_score_string << child->exafs_score << ",";
			potential_energy_string << child->potential_energy << ",";
		}
	}

	this->output_stream << exafs_score_string.str() << std::endl;
	this->output_stream << potential_energy_string.str() << std::endl;
}

void EXAFSMOODE::finalStats() {

	std::cout << "Final Stats." << std::endl;

	this->output_stream.close();

	std::stringstream output_chromosome_filename;
	int counter = 0;
	for (std::vector<Chromosome>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		std::cout << "Rank == " << child->rank << std::endl;
		if (child->rank == 1) {
			
			output_chromosome_filename << "/best-chromosome-" << (counter++) << ".pdb";
			std::cout << "Writing out: " << this->stats_folder << output_chromosome_filename.str() << std::endl;
			
			this->exafs_evaluator->updateAtoms(child->atoms);
			this->exafs_evaluator->writePDB(this->stats_folder + output_chromosome_filename.str());
			output_chromosome_filename.str("");
		}
	}
	
}

bool EXAFSMOODE::convergence() {

	std::vector<Chromosome> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort_moode);
	std::vector<Chromosome>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) == 1;
}

bool EXAFSMOODE::convergence(double rate) {

	std::vector<Chromosome> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort_moode);
	std::vector<Chromosome>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) <= (rate * this->population.size());
}

void EXAFSMOODE::rankPopulation() {

	int currentRank = 1;

	std::vector<Chromosome> current_population = this->population;
	std::vector<Chromosome> new_population;
	// Reset all ranks.
	for (int i = 0; i < (int)current_population.size(); ++i) {
		current_population[i].rank = -1;
	}

	while((int)current_population.size() > 0) {

		int currentRankCount = 0;
		for (int i = 0; i < (int)current_population.size(); ++i) {
			
			bool should_keep = true;
			for (int j = 0; j < (int)current_population.size(); ++j) {
				
				if (i == j) continue;

				int result = this->chromosomeDominates(current_population[i],current_population[j]);

				if (result == -1) {
				
					if (this->chromosomeDominates(current_population[j],current_population[i]) == 1) {
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

		std::vector<Chromosome> unranked_population;

		for (int i = 0; i < (int)current_population.size(); ++i) {
			if (current_population[i].rank != -1) {
				new_population.push_back(current_population[i]);
			} else {
				unranked_population.push_back(current_population[i]);
			}
		}

		current_population = unranked_population;
	}

	std::cout << "Population size = " << new_population.size() << std::endl;
	this->population = new_population;
	
}

/*
	Checks is Chromosome a dominate Chromosome b

	Returns:

	-1: no
	0: equal
	1: yes
*/
int EXAFSMOODE::chromosomeDominates(Chromosome& a, Chromosome& b) {

	bool at_least_one_better = false;

	if (a.exafs_score <= b.exafs_score) {
		if (a.exafs_score < b.exafs_score) {
			at_least_one_better = true;
		}
	} else {
		return -1;
	}

	if (a.potential_energy <= b.potential_energy) {
		if (a.potential_energy < b.potential_energy) {
			at_least_one_better = true;
		}
	} else {
		return -1;
	}

	return at_least_one_better ? 1 : 0;
}