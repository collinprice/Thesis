#include <string>
#include <vector>
#include "atom.h"

class Molecule {

	public:
		std::vector<Atom> atoms;
		double cost;
		bool isEvaluated;

		Molecule(std::vector<Atom> atoms);
		// Molecule operator= (Molecule);
		int size();

		// Utility
		bool write_to_disk(std::string filename);
		static Molecule read_from_disk(std::string filename);
};