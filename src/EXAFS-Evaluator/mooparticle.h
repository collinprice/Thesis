#include "chromosome.h"
#include "particlearchive.h"
// #include "basicparticle.h"

// class ParticleArchive;

class MOOParticle : public Chromosome {

	public:
		ParticleArchive archive;

		std::vector< PDBAtom > velocity;

		MOOParticle();
		MOOParticle(double range, std::vector<PDBAtom> atoms);
		MOOParticle(const MOOParticle& other);
		// ~MOOParticle();

		void updateBest();
		void updatePosition();
		void updateVelocity(BasicParticle global_best, double inertia, double social, double cognitive);

		BasicParticle getBasic();
	
	private:
		double unifRand();
		double randomRange(double max);
};