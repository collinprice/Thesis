#include "exafsevaluator.h"
#include "mooparticle.h"
// #include "basicparticle.h"
// #include "particlearchive.h"

#include <vector>
#include <fstream>

class MOOPSO {
	public:

		MOOPSO(EXAFSEvaluator* exafs_evaluator, double inertia, double social, double cognitive, double velocity_range, int max_generations, std::string results_file);
		~MOOPSO();
		void begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations);

	private:

		std::vector<MOOParticle> population;

		ParticleArchive archive;
		// MOOParticle global_best_particle;
		BasicParticle global_best_particle;

		EXAFSEvaluator* exafs_evaluator;
		double inertia;
		double social;
		double cognitive;
		double velocity_range;
		
		int max_generations;
		std::string results_file;
		std::ofstream output_stream;
		std::string stats_folder;

		void initPopulation(std::vector< std::vector<PDBAtom> > population);

		void updateVelocities();
		void updatePositions();
		void updateGlobalBest();

		void evaluate( MOOParticle& child );
		void evaluatePopulation();

		// Particle best_particle();
		// void saveBestParticles();

		double unifRand();
		double randInt(int max);

		void initStats();
		void recordStats();
		void finalStats();

		bool convergence();
		bool convergence(double rate);
};